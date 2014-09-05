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

#ifndef TESTMONODOMAINGENERATEREFERENCESOLUTIONS_HPP_
#define TESTMONODOMAINGENERATEREFERENCESOLUTIONS_HPP_

/*
 * = Create accurate reference monodomain tissue simulation solutions =
 *
 * On this wiki page we describe in detail the code that is used to run this example from the paper.
 *
 * We run a 1D monodomain simulation on a highly refined PDE time step with a handful of the models.
 *
 * == Code overview ==
 *
 * The first thing to do is to include the necessary header files.
 */

#include <cxxtest/TestSuite.h>

#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>

// Headers specific to this project
#include "CellModelUtilities.hpp"
#include "DynamicModelCellFactory.hpp"
#include "MonodomainProblem.hpp"

// This header is needed to allow us to run in parallel
#include "PetscSetupAndFinalize.hpp"

class TestMonodomainGenerateReferenceSolutionsLiteratePaper : public CxxTest::TestSuite
{
public:
    void TestMonodomainGenerateReferenceSolutions() throw (Exception)
    {
        std::vector<FileFinder> models = CellModelUtilities::GetListOfModels();

        std::vector<std::string> pycml_options;
        pycml_options.push_back("--cvode");
        pycml_options.push_back("--expose-annotated-variables");

        DynamicModelCellFactory cell_factory(models[0], pycml_options);

        /* We will auto-generate a mesh this time, and pass it in, rather than
         * provide a mesh file name. This is how to generate a cuboid mesh with
         * a given spatial stepsize h.
         *
         * Using a `DistributedTetrahedralMesh` is faster than `TetrahedralMesh` when running on multiple processes.
         * However, it permutes the node ordering for output. Most of time time this won't matter, but later in this
         * test we want to access specific node indices. One method of doing this is to ask `HeartConfig` to use the
         * original node ordering for the output.
         *
         */
        DistributedTetrahedralMesh<1,1> mesh;
        double h=0.01;
        mesh.ConstructRegularSlabMesh(h, 1 /*length*/);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);

        /*
         * EMPTYLINE
         *
         * Set the simulation duration, etc, and create an instance of the cell factory.
         * One thing that should be noted for monodomain problems, the ''intracellular
         * conductivity'' is used as the monodomain effective conductivity (not a
         * harmonic mean of intra and extracellular conductivities). So if you want to
         * alter the monodomain conductivity call
         * `HeartConfig::Instance()->SetIntracellularConductivities`
         */
        HeartConfig::Instance()->SetSimulationDuration(1000); //ms
        HeartConfig::Instance()->SetOutputDirectory("Frontiers/MonodomainReference/" + models[0].GetLeafNameNoExtension() + "/results");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);

        /* Now we declare the problem class */
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        /* If a mesh-file-name hasn't been set using `HeartConfig`, we have to pass in
         * a mesh using the `SetMesh` method (must be called before `Initialise`). */
        monodomain_problem.SetMesh(&mesh);

        /* `SetWriteInfo` is a useful method that means that the min/max voltage is
         * printed as the simulation runs (useful for verifying that cells are stimulated
         * and the wave propagating, for example) (although note scons does buffer output
         * before printing to screen) */
        monodomain_problem.SetWriteInfo();

        /* Finally, call `Initialise` and `Solve` as before */
        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        /* This part is just to check nothing has accidentally been changed in this example */
        ReplicatableVector voltage(monodomain_problem.GetSolution());
        TS_ASSERT_DELTA(voltage[0], 34.9032, 1e-2);
    }
};

#endif // TESTMONODOMAINGENERATEREFERENCESOLUTIONS_HPP_
