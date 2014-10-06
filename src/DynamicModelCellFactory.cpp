/*

Copyright (c) 2005-2014, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "DynamicModelCellFactory.hpp"

#include <sstream>
#include "HeartFileFinder.hpp"
#include "CellMLToSharedLibraryConverter.hpp"
#include "AbstractCardiacCellInterface.hpp"
#include "Warnings.hpp"
// This is needed to prevent the chaste_libs=0 build failing
// on tests that use a dynamically loaded CVODE model
#include "AbstractCvodeCell.hpp"

DynamicModelCellFactory::DynamicModelCellFactory(const FileFinder& rModelFile,
                                                 OutputFileHandler& rOutputDir,
                                                 Solvers::Value solver,
                                                 bool useLookupTables,
                                                 bool highTol)
    : AbstractCardiacCellFactory<1u>(),
      mStrictTolerance(highTol),
      mpLookupTables(NULL),
      mSolver(solver)
{
    // Create a simple stimulus to use in the CreateCardiacCellForTissueNode method.
    // The stimulus is a bit of a compromise, needs to be big enough to fire off all the models, but not
    // so big that it pushes the voltage up to 'unphysiological' for certain models.
    mpSimpleStimulus = boost::shared_ptr<SimpleStimulus>(new SimpleStimulus(-30000, 4, 1)); // magnitude, duration, start time.

    // Do a single run to do the conversion and create the .so model that we will want to use (discard model it gives).
    CellModelUtilities::CreateCellModel(rModelFile,rOutputDir,solver,useLookupTables);

    mSharedLibraryLocation = rOutputDir.FindFile("lib" + rModelFile.GetLeafNameNoExtension() + ".so");
    assert(mSharedLibraryLocation.IsFile());
}

DynamicModelCellFactory::~DynamicModelCellFactory()
{
}

AbstractCardiacCellInterface* DynamicModelCellFactory::CreateCardiacCellForTissueNode(Node<1u>* pNode)
{
    AbstractCardiacCellInterface* p_cell = NULL;

#ifndef CHASTE_CAN_CHECKPOINT_DLLS
    if (HeartConfig::Instance()->GetCheckpointSimulation())
    {
        EXCEPTION("Checkpointing is not compatible with dynamically loaded cell models on Boost<1.37.");
    }
#endif // CHASTE_CAN_CHECKPOINT_DLLS

    //
    // Load model from shared library
    //

    CellMLToSharedLibraryConverter converter;
    DynamicCellModelLoaderPtr p_loader = converter.Convert(mSharedLibraryLocation, false);

    // Load up the cell with a default solver
    p_cell = p_loader->CreateCell(this->mpSolver, this->mpZeroStimulus);

    // Overwrite the solver
    CellModelUtilities::SetCellModelSolver(p_cell, mSolver);

    /*
     * If this is anywhere inside x =< 0.1cm we apply a stimulus, otherwise we don't.
     *
     * N.B. Just applying a stimulus at one node on the edge is a bit numerically fragile
     * as this looks quite different to the next node on different mesh spacings.
     */
    if (pNode->rGetLocation()[0] <= 0.1)
    {
        //std::cout << "Generating cell model with a stimulus" << std::endl << std::flush;
        p_cell->SetStimulusFunction(this->mpSimpleStimulus);
    }

    // Generate lookup tables if present
    mpLookupTables = p_cell->GetLookupTableCollection();

    if (mStrictTolerance && dynamic_cast<AbstractCvodeCell*>(p_cell))
    {
        static_cast<AbstractCvodeCell*>(p_cell)->SetTolerances(1e-7,1e-9);
        // NB. We can't set the CVODE cells to 'no reset' here, as it would be overridden by the
        // Abstract cell factory...!
    }

    return p_cell;
}

void DynamicModelCellFactory::FreeLookupTableMemory()
{
    if (mpLookupTables)
    {
        mpLookupTables->FreeMemory();
    }
}
