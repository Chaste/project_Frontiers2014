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
 * = Create reference monodomain tissue simulation solutions with CVODE =
 *
 * On this wiki page we describe in detail the code that is used to run this example from the paper.
 *
 * We run a 1D monodomain simulation on a suitably refined PDE time step with a handful of the models.
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

// Chaste 'heart' headers
#include "MonodomainProblem.hpp"
#include "CellProperties.hpp"
#include "AbstractCvodeCell.hpp"

// This header is needed to allow us to run in parallel
#include "PetscSetupAndFinalize.hpp"

class TestMonodomainGenerateReferenceSolutionsLiteratePaper : public CxxTest::TestSuite
{
public:
    void TestMonodomainGenerateReferenceSolutions() throw (Exception)
    {
        // We don't want to find this confusing matters!!
        EXIT_IF_PARALLEL;

        /*
         * This test was run with the following values inserted here
         *  * 0.01 ms (typically a fine timestep)
         *  * 0.1 ms (about average for most studies)
         */
        std::vector<double> pde_timesteps = boost::assign::list_of(0.1)(0.01);

        /*
         * This test is run with the following values, selected by looking at the output of
         * the convergence test.
         */
        double h = 0.01;

        std::vector<FileFinder> all_models = CellModelUtilities::GetListOfModels();

        // A list of models that we want to do tissue simulations with.
        std::vector<std::string> models_to_use = boost::assign::list_of("luo_rudy_1991")
                                                 ("beeler_reuter_model_1977")
                                                 ("nygren_atrial_model_1998")
                                                 ("ten_tusscher_model_2004_epi")
                                                 ("grandi_pasqualini_bers_2010_ss")
                                                 ("shannon_wang_puglisi_weber_bers_2004")
                                                 ("iyer_model_2007");

        OutputFileHandler overall_results_folder("Frontiers/MonodomainReference/", true); // Wipe!
        out_stream p_file = overall_results_folder.OpenOutputFile("error_summary_tissue.txt");
        *p_file << std::setiosflags(std::ios::scientific) << std::setprecision(8);

        /* Repository data location */
        FileFinder this_file(__FILE__);
        FileFinder repo_data("data", this_file);

        // Loop over models
        BOOST_FOREACH(std::string model, models_to_use)
        {
            // Find the FileFinder associated with the model we want.
            FileFinder model_to_use;
            for (unsigned i=0; i<all_models.size(); i++)
            {
                if (all_models[i].GetLeafNameNoExtension()==model)
                {
                    model_to_use = all_models[i];
                    break;
                }
            }

            std::stringstream output_folder_stream;
            output_folder_stream << "Frontiers/MonodomainReference/" << model;
            std::string output_folder = output_folder_stream.str();

            // Use a different constructor
            OutputFileHandler base_model_handler(output_folder);
            DynamicModelCellFactory cell_factory(model_to_use,
                                                 base_model_handler,
                                                 Solvers::CVODE_ANALYTIC_J,
                                                 false,// Whether to use lookup tables.
                                                 false); // false for using CVODE with regular tolerances

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
             * `HeartConfig::Instance()->SetIntracellularConductivities()`
             */
            HeartConfig::Instance()->SetSimulationDuration(500); //ms
            HeartConfig::Instance()->SetOutputDirectory(output_folder + "/results");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");
            HeartConfig::Instance()->SetVisualizeWithVtk(true);

            BOOST_FOREACH(double pde_timestep, pde_timesteps)
            {

                HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(pde_timestep, pde_timestep, 0.1);

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

                HeartEventHandler::Reset();

                /* Finally, call `Initialise` and `Solve` as usual */
                monodomain_problem.Initialise();

                try
                {
                    monodomain_problem.Solve();
                }
                catch (Exception& e)
                {
                    WARNING("Model '" << model << "' simulation failed with: " << e.GetMessage());
                    continue;
                }
                HeartEventHandler::Headings();
                HeartEventHandler::Report();

                /*
                 * Read some of the outputted data back in, and evaluate AP properties at the last node,
                 * as per the single cell stuff.
                 */
                Hdf5DataReader data_reader = monodomain_problem.GetDataReader();
                std::vector<double> times = data_reader.GetUnlimitedDimensionValues();

                std::vector<double> last_node = data_reader.GetVariableOverTime("V", mesh.GetNumNodes()-1u);

                try
                {
                    std::vector<double> errors = CellModelUtilities::GetTissueErrors(times, last_node, model, pde_timestep);

                    std::cout << "Model: " << model << " Square error = " << errors[0] << ", APD90 error = " << errors[1] <<
                            ", APD50 error = " << errors[2] << ", APD30 error = " << errors[3] << ", V_max error = " << errors[4] <<
                            ", V_min error = " << errors[5] << ", dVdt_max error = " << errors[6] << std::endl;

                    // Write to file too.
                    *p_file << model << "\t" << pde_timestep;
                    for (unsigned i=0; i<errors.size(); i++)
                    {
                        *p_file << "\t" << errors[i];
                    }
                    *p_file << std::endl;
                }
                catch(Exception &e)
                {
                    WARNING("Model " << model << ": analysis of voltage at last node failed.");
                    continue;
                }
            }
        }

        p_file->close();

        // Now copy it into the repository.
        FileFinder ref_data = overall_results_folder.FindFile("error_summary_tissue.txt");
        ref_data.CopyTo(repo_data);
    }
};

#endif // TESTMONODOMAINGENERATEREFERENCESOLUTIONS_HPP_
