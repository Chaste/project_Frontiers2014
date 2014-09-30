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

#ifndef TESTMONODOMAINCALCULATEREQUIREDTIMESTEPS_HPP_
#define TESTMONODOMAINCALCULATEREQUIREDTIMESTEPS_HPP_

/*
 * = Calculate the required ODE solver timesteps to meet target PDE accuracy =
 *
 * On this wiki page we describe in detail the code that is used to run this example from the paper.
 *
 * This test performs a monodomain simulation with a range of models, using with each ODE solver.
 * It reduces the time step used until the error metric (based on sum of squares for voltage trace at last node)
 * is within tolerances based on the reference solutions from
 * in PaperTutorials/Frontiers2014/MonodomainConvergence
 *
 * This is intended to define a time step required for each method to
 * get a numerical solution of comparable accuracy, for fair timing comparisons.
 *
 * == Code overview ==
 *
 * The first thing to do is to include the necessary header files.
 */

// The testing framework we use
#include <cxxtest/TestSuite.h>

// Standard C++ libraries
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#include <boost/pointer_cast.hpp>

// Headers specific to this project
#include "CellModelUtilities.hpp"
#include "DynamicModelCellFactory.hpp"

// Chaste 'heart' headers
#include "MonodomainProblem.hpp"
#include "CellProperties.hpp"
#include "AbstractCvodeCell.hpp"

// This header is needed to allow us to run in parallel
#include "PetscSetupAndFinalize.hpp"

class TestMonodomainCalculateRequiredTimestepsLiteratePaper : public CxxTest::TestSuite
{
public:
    void TestMonodomainCalculateTimesteps() throw (Exception)
    {
        // We don't want to find this confusing matters!!
        EXIT_IF_PARALLEL;

        const double required_mrms_error = 0.05; // 5%

        /*
         * This test was run with the following values inserted here
         *  * 0.01 ms (typically a fine timestep)
         *  * 0.1 ms (about average for most studies)
         */
        std::vector<double> pde_timesteps = boost::assign::list_of(0.1)(0.01);

        /*
         * This test is run with the following values, selected by looking at the output of
         * [TestMonodomainConvergenceLiteratePaper the convergence test].
         */
        double h = 0.01;

        std::vector<FileFinder> all_models = CellModelUtilities::GetListOfModels();

        /* A list of models that we want to do tissue simulations with, and then the solvers to test. */
        std::vector<std::string> models_to_use = boost::assign::list_of("luo_rudy_1991")
                                                 ("beeler_reuter_model_1977")
                                                 ("nygren_atrial_model_1998")
                                                 ("ten_tusscher_model_2004_epi")
                                                 ("grandi_pasqualini_bers_2010_ss")
                                                 ("shannon_wang_puglisi_weber_bers_2004")
                                                 ("iyer_model_2007");

        std::vector<Solvers::Value> solvers = boost::assign::list_of
                (Solvers::CVODE_ANALYTIC_J)(Solvers::CVODE_NUMERICAL_J)
                (Solvers::FORWARD_EULER)(Solvers::BACKWARD_EULER)
                (Solvers::RUNGE_KUTTA_2)(Solvers::RUNGE_KUTTA_4)
                (Solvers::RUSH_LARSEN)(Solvers::GENERALISED_RUSH_LARSEN_1)(Solvers::GENERALISED_RUSH_LARSEN_2);

        OutputFileHandler overall_results_folder("Frontiers/MonodomainCalculateTimesteps/", false); // Don't wipe!
        out_stream p_file = overall_results_folder.OpenOutputFile("required_steps_tissue.txt");
        *p_file << std::setiosflags(std::ios::scientific) << std::setprecision(8);

        /* Repository data location */
        FileFinder this_file(__FILE__);
        FileFinder repo_data("data", this_file);
        FileFinder model_folder("../cellml", this_file);

        // Loop over models
        BOOST_FOREACH(std::string model, models_to_use)
        {
            /* Find the CellML file for this model. */
            FileFinder model_to_use(model + ".cellml", model_folder);

            /* Iterate over each available solver, using a handy boost method. */
            BOOST_FOREACH(Solvers::Value solver, solvers)
            {
                bool cvode_solver = ((solver==Solvers::CVODE_ANALYTIC_J) || (solver==Solvers::CVODE_NUMERICAL_J));

                std::cout << "Calculating timestep for " << model << " with solver " << CellModelUtilities::GetSolverName(solver) << std::endl;

                std::stringstream output_folder_stream;
                output_folder_stream << "Frontiers/MonodomainCalculateTimesteps/" << model << "/" << solver;
                std::string output_folder = output_folder_stream.str();

                // Use a different constructor
                OutputFileHandler base_model_handler(output_folder);
                DynamicModelCellFactory cell_factory(model_to_use,
                                                     base_model_handler,
                                                     solver,
                                                     false); // Whether to use lookup tables.

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

                BOOST_FOREACH(double pde_timestep, pde_timesteps)
                {
                    /*
                     * Set the simulation duration, etc, and create an instance of the cell factory.
                     * One thing that should be noted for monodomain problems, the ''intracellular
                     * conductivity'' is used as the monodomain effective conductivity (not a
                     * harmonic mean of intra and extracellular conductivities). So if you want to
                     * alter the monodomain conductivity call
                     * `HeartConfig::Instance()->SetIntracellularConductivities()`
                     */
                    HeartConfig::Instance()->SetSimulationDuration(500); //ms

                    unsigned timestep_divisor = 1u;
                    unsigned refinement_idx = 1u;
                    bool within_tolerance = false;
                    double sampling_time = 0.1; //ms

                    /* Keep solving with smaller timesteps until we meet the desired accuracy */
                    double elapsed_time;
                    do
                    {
                        double ode_timestep = pde_timestep / timestep_divisor;
                        std::stringstream refinement_description_stream;
                        if (cvode_solver)
                        {
                            // Reset time step to be equal to the PDE timestep
                            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(pde_timestep, pde_timestep, sampling_time);
                            if (refinement_idx>=7u)
                            {
                                // CVODE is struggling, probably because it is hitting a singularity?
                                break;
                            }
                            ode_timestep = (double)(refinement_idx); // We just hijack this variable for the logs.
                            refinement_description_stream << " CVODE tolerances of " << pow(10,-1.0-refinement_idx) <<
                                    ", " << pow(10,-3.0-refinement_idx) ;
                        }
                        else
                        {
                            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(ode_timestep, pde_timestep, sampling_time);
                            refinement_description_stream << " timestep of " << ode_timestep << "ms";
                        }
                        std::string refinement_description = refinement_description_stream.str();

                        std::cout << std::endl << "Computing error for PDE dt = " << pde_timestep << " model: " << model << " solver: '" << CellModelUtilities::GetSolverName(solver) << "' with" << refinement_description << std::endl;

                        std::stringstream output_subfolder;
                        output_subfolder << "/results_pde_" <<  pde_timestep << "_ode_" << ode_timestep;

                        HeartConfig::Instance()->SetOutputDirectory(output_folder + output_subfolder.str());
                        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
                        HeartConfig::Instance()->SetVisualizeWithVtk(true);

                        /* Now we declare the problem class */
                        MonodomainProblem<1> monodomain_problem( &cell_factory );

                        /* If a mesh-file-name hasn't been set using `HeartConfig`, we have to pass in
                         * a mesh using the `SetMesh` method (must be called before `Initialise`). */
                        monodomain_problem.SetMesh(&mesh);

                        HeartEventHandler::Reset();

                        /* Finally, call `Initialise` and `Solve` as usual */
                        monodomain_problem.Initialise();

                        /* The cells should now be set up, we have to 'hack in' and alter their tolerances */
                        if (cvode_solver)
                        {
                            /* Not that we have special methods to iterate over nodes of a distributed mesh
                             * which will work in parallel settings too. */
                            DistributedVectorFactory* p_factory = mesh.GetDistributedVectorFactory();
                            Vec monodomain_vec = p_factory->CreateVec();
                            DistributedVector monodomain_distributed_vec = p_factory->CreateDistributedVector(monodomain_vec);

                            for (DistributedVector::Iterator index=monodomain_distributed_vec.Begin();
                                 index != monodomain_distributed_vec.End();
                                 ++index)
                            {
                                AbstractCardiacCellInterface* p_cell = monodomain_problem.GetMonodomainTissue()->GetCardiacCell(index.Global);
                                CellModelUtilities::SetCvodeTolerances(p_cell, refinement_idx);
                            }
                        }

                        timestep_divisor *= 2; // Ready for next iteration
                        refinement_idx++; // Ready for next iteration

                        try
                        {
                            Timer::Reset();
                            monodomain_problem.Solve();
                            elapsed_time = Timer::GetElapsedTime();
                        }
                        catch (Exception& e)
                        {
                            std::cout << model << " failed to solve with ODE timestep/refinement of " << ode_timestep << ", got: " << e.GetMessage() << std::endl << std::flush;
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

                        /* Put this trace into a file for easy plotting and comparison with the reference later on. */
                        OutputFileHandler handler(output_folder + output_subfolder.str(), false);
                        out_stream p_trace_file = handler.OpenOutputFile("last_node_trace.dat");
                        for (unsigned i=0; i<times.size(); i++)
                        {
                            *p_trace_file << times[i] << "\t" << last_node[i] << std::endl;
                        }
                        p_trace_file->close();

                        /* Analyse the error associated with this run compared to the reference trace in the repository. */
                        try
                        {
                            std::vector<double> errors = CellModelUtilities::GetTissueErrors(times, last_node, model, pde_timestep);

                            std::cout << "Model: " << model << " Square error = " << errors[0] << ", APD90 error = " << errors[1] <<
                                    ", APD50 error = " << errors[2] << ", APD30 error = " << errors[3] << ", V_max error = " << errors[4] <<
                                    ", V_min error = " << errors[5] << ", dVdt_max error = " << errors[6] << ", MRMS error = " << errors[7] << std::endl;

                            within_tolerance = (errors[7] <= required_mrms_error);

                            std::cout << "For timestep " << ode_timestep << " ms MRMS error = " << errors[7] << ". Good enough ? " << within_tolerance << std::endl << std::flush;

                            // Write to file too.
                            *p_file << model << "\t" << pde_timestep << "\t" << solver << "\t" << ode_timestep;
                            for (unsigned i=0; i<errors.size(); i++)
                            {
                                *p_file << "\t" << errors[i];
                            }
                            *p_file << "\t" << within_tolerance << "\t" << elapsed_time;
                            *p_file << std::endl;
                        }
                        catch (const Exception& r_e)
                        {
                            std::cout << "Exception computing error for model " << model << " solver " << CellModelUtilities::GetSolverName(solver);
                            std::cout << r_e.GetMessage() << std::endl;
                            WARNING(r_e.GetMessage());
                        }

                        /* We stop bothering if the simulation took longer than 5 minutes
                         * (as we know even Iyer on PDE 0.01ms was solved in 140 seconds using CVODE)
                         * so we're in a regime where it is less accurate and takes at least twice as long
                         * as CVODE, we don't need to know any more to discount this! */
                        if (elapsed_time >= 300)
                        {
                            std::stringstream message;
                            message << "Model: " << model << " with solver '" << CellModelUtilities::GetSolverName(solver) << "' and timestep " << ode_timestep << " took longer than 5 minutes to solve, so we're not refining any more.\n";
                            std::cout << message.str() << std::flush;
                            WARNING(message.str());
                            break;
                        }
                    }
                    while (!within_tolerance && timestep_divisor <= 2048);
                    /* The above means we allow at most 12 refinements of the time step (2^12^ = 2048).*/
                }

                // Free memory for lookup tables if used
                cell_factory.FreeLookupTableMemory();
            }
        }

        p_file->close();

        /* Copy to repository for storage and use by the timing tests */
        FileFinder ref_data = overall_results_folder.FindFile("required_steps_tissue.txt");
        ref_data.CopyTo(repo_data);
    }
};

#endif // TESTMONODOMAINCALCULATEREQUIREDTIMESTEPS_HPP_
