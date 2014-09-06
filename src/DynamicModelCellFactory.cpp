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
                                                 const std::vector<std::string>& rPyCmlOptions)
    : AbstractCardiacCellFactory<1u>()
{
    // Create a simple stimulus to use.
    mpSimpleStimulus = boost::shared_ptr<SimpleStimulus>(new SimpleStimulus(-250000, 3, 1)); // magnitude, duration, start time.

    // Make a new output folder for this model.
    mpHandler = new OutputFileHandler("Frontiers/MonodomainReference/" + rModelFile.GetLeafNameNoExtension(), true);

    // First let's copy the CellML file to this folder to make the
    // .so and outputted files end up there too.
    FileFinder copied_cellml_file = mpHandler->CopyFileTo(rModelFile);

    // Also copy the .out file across, if it exists.
    FileFinder maple_output_file(rModelFile.GetLeafNameNoExtension() + ".out", rModelFile);
    if (maple_output_file.IsFile())
    {
        mpHandler->CopyFileTo(maple_output_file);
    }

    // We first collectively convert this CellML file to shared library,
    // which each process can load to assign a new model to its nodes.
    DynamicCellModelLoaderPtr p_loader = LoadDynamicModel(copied_cellml_file, rPyCmlOptions, true);

    mSharedLibraryLocation = FileFinder(p_loader->GetLoadableModulePath(), RelativeTo::Absolute);
    std::cout << "Shared library location is " << mSharedLibraryLocation.GetAbsolutePath() << "\n";
    assert(mSharedLibraryLocation.Exists());
}

DynamicCellModelLoaderPtr DynamicModelCellFactory::LoadDynamicModel(
        const FileFinder& rModelFile,
        const std::vector<std::string>& rPyCmlOptions,
        bool isCollective)
{
    assert(rModelFile.Exists());
    CellMLToSharedLibraryConverter converter(true); // true to preserve generated source for inspection

    /* We only want to write out the options file once at the beginning */
    if (isCollective)
    {
        converter.CreateOptionsFile(*mpHandler,
                                    rModelFile.GetLeafNameNoExtension(),
                                    rPyCmlOptions);
    }
    return converter.Convert(rModelFile, isCollective);
}

DynamicModelCellFactory::~DynamicModelCellFactory()
{
    delete mpHandler;
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

    // Load model from shared library

    // options won't be used with a last 'false' argument, so we pass in an empty vec.
    std::vector<std::string> empty_options;

    DynamicCellModelLoaderPtr p_loader = LoadDynamicModel(mSharedLibraryLocation, empty_options, false);

    /* If this is the first node at x=0 we apply a stimulus, otherwise we don't */
    if (pNode->rGetLocation()[0] < 1e-6)
    {
        std::cout << "Generating cell model with a stimulus\n";
        p_cell = p_loader->CreateCell(this->mpSolver, this->mpSimpleStimulus);
    }
    else
    {
        p_cell = p_loader->CreateCell(this->mpSolver, this->mpZeroStimulus);
    }

    // Generate lookup tables if present
    p_cell->GetLookupTableCollection();
    return p_cell;
}


