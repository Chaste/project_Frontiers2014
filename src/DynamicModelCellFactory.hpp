/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef DYNAMICMODELCELLFACTORY_HPP_
#define DYNAMICMODELCELLFACTORY_HPP_

#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

#include "AbstractCardiacCellFactory.hpp"
#include "CellModelUtilities.hpp"
#include "SimpleStimulus.hpp"

#include "DynamicCellModelLoader.hpp"

/**
 * This is a cardiac cell factory which is a stripped down version of
 * HeartConfigRelatedCellFactory which has the ability to dynamically
 * convert CellML files into C++ code, compile this, and load the resulting model.
 */
class DynamicModelCellFactory : public AbstractCardiacCellFactory<1u>
{
private:
    /** Whether to switch CVODE onto a very high tolerance for the reference solution. */
    bool mStrictTolerance;

    /** A simple stimulus to apply at one end at t=1.0 */
    boost::shared_ptr<SimpleStimulus> mpSimpleStimulus;

    /** The shared library location. */
    FileFinder mSharedLibraryLocation;

    /** The lookup tables used by cells we create. */
    AbstractLookupTableCollection* mpLookupTables;

    /** The solver type we are using */
    Solvers::Value mSolver;

public:
    /**
     * A constructor which will use CellModelUtilities to dynamically convert a CellML
     * file into C++ code, compile it and load the resulting model object.
     *
     * @param rModelFile  the CellML file to convert
     * @param rOutputDir  handler for the folder in which to create generated code
     * @param solver  which cell model solver to use
     * @param useLookupTables  whether to use lookup tables to speed up simulations
     * @param highTol  Whether we are generating the reference solution.
     */
    DynamicModelCellFactory(const FileFinder& rModelFile,
                            OutputFileHandler& rOutputDir,
                            Solvers::Value solver,
                            bool useLookupTables,
                            bool highTol=false);

    /** Destructor*/
    ~DynamicModelCellFactory();

    /**
     * @return a newly created correct stimulated tissue cell for a given region in the mesh
     * @param pNode pointer to the node.
     */
    AbstractCardiacCellInterface* CreateCardiacCellForTissueNode(Node<1u>* pNode);

    /**
     * This can be called at the end of a simulation to free the memory used by lookup tables, if any.
     */
    void FreeLookupTableMemory();
};


#endif /*DYNAMICMODELCELLFACTORY_HPP_*/
