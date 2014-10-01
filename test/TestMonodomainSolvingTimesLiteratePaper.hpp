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

#ifndef TESTMONODOMAINSOLVINGTIMESLITERATEPAPER_HPP_
#define TESTMONODOMAINSOLVINGTIMESLITERATEPAPER_HPP_

/*
 * = Benchmark monodomain tissue simulation solving times with different numerical methods =
 *
 * On this wiki page we describe in detail the code that is used to run this example from the paper.
 *
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

// Chaste headers
#include "MonodomainProblem.hpp"
#include "CellProperties.hpp"
#include "Timer.hpp"
#include "ChasteBuildRoot.hpp"

// This header is needed to allow us to run in parallel
#include "PetscSetupAndFinalize.hpp"

class TestMonodomainSolvingTimesLiteratePaper : public CxxTest::TestSuite
{
public:
    void TestMonodomainSolvingTimes() throw (Exception)
    {
        /*
         * Load the archived timesteps that are appropriate for each model/solver/pde_timestep combo out of the
         * file test/data/required_timesteps_tissue.txt
         *
         * We load and print them to screen in a small separate test here, so it can be run easily
         */
        LoadTimestepFile();

        /*
         * This test was run with the following values inserted here
         *  * 0.01 ms (typically a fine timestep)
         *  * 0.1 ms (about average for most studies)
         */
        std::vector<double> pde_timesteps = boost::assign::list_of(0.1)(0.01);

        /*
         * This test is run with the following values, selected by looking at the output of
         * [MonodomainConvergence the convergence test].
         */
        double h = 0.01;

        // A list of models that we want to do tissue simulations with.
        std::vector<std::string> models_to_use = boost::assign::list_of("luo_rudy_1991")
                                                 ("beeler_reuter_model_1977")
                                                 ("nygren_atrial_model_1998")
                                                 ("ten_tusscher_model_2004_epi")
                                                 ("grandi_pasqualini_bers_2010_ss")
                                                 ("shannon_wang_puglisi_weber_bers_2004")
                                                 ("iyer_model_2007");

        std::vector<Solvers::Value> solvers = boost::assign::list_of(Solvers::CVODE_ANALYTIC_J)(Solvers::CVODE_NUMERICAL_J)
                (Solvers::FORWARD_EULER)(Solvers::BACKWARD_EULER)
                (Solvers::RUNGE_KUTTA_2)(Solvers::RUNGE_KUTTA_4)
                (Solvers::RUSH_LARSEN)(Solvers::GENERALISED_RUSH_LARSEN_1)(Solvers::GENERALISED_RUSH_LARSEN_2);

        OutputFileHandler overall_results_folder("Frontiers/MonodomainTimings/", true); // Wipe!
        out_stream p_file;
        if (PetscTools::AmMaster())
        {
            p_file = overall_results_folder.OpenOutputFile(ChasteBuildType() + "_timings_tissue.txt");
            *p_file << std::setiosflags(std::ios::scientific) << std::setprecision(8);
        }

        /* Repository data location */
        FileFinder this_file(__FILE__);
        FileFinder repo_data("data", this_file);
        FileFinder model_folder("../cellml", this_file);

        // Loop over models
        BOOST_FOREACH(std::string model, models_to_use)
        {
            /* Find the CellML file for this model. */
            FileFinder model_to_use(model + ".cellml", model_folder);

            /* Iterate over each available solver, using a handy boost method */
            BOOST_FOREACH(Solvers::Value solver, solvers)
            {
                // Used to get the correct timestep to use later on...
                std::pair<std::string, Solvers::Value> model_solver_combo(model,solver);

                std::vector<bool> lookup_table_options = boost::assign::list_of(false)(true);
                BOOST_FOREACH(bool use_lookup_tables, lookup_table_options)
                {
                    std::string lookup_tables_description = "";
                    if (use_lookup_tables)
                    {
                        lookup_tables_description = " with lookup tables";
                    }
                    std::cout << "Running timings for " << model << " with solver " << CellModelUtilities::GetSolverName(solver) << lookup_tables_description << std::endl;

                    std::stringstream output_folder_stream;
                    output_folder_stream << "Frontiers/MonodomainTimings/" << model << "/" << CellModelUtilities::GetSolverName(solver) << "_" << use_lookup_tables;
                    std::string output_folder = output_folder_stream.str();

                    // Use a different constructor
                    OutputFileHandler base_model_handler(output_folder);
                    DynamicModelCellFactory cell_factory(model_to_use,
                                                         base_model_handler,
                                                         solver,
                                                         use_lookup_tables);

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

                        double ode_timestep = pde_timestep; // Default for CVODE solvers
                        if ((solver != Solvers::CVODE_ANALYTIC_J) && (solver != Solvers::CVODE_NUMERICAL_J))
                        {
                            bool found = false;
                            bool accurate_enough = false;
                            if (pde_timestep==0.1)
                            {
                                std::map<std::pair<std::string, Solvers::Value>, std::pair<double,bool> >::iterator it = mTimestepsPde0_1.find(model_solver_combo);
                                if (it != mTimestepsPde0_1.end())
                                {
                                    ode_timestep = (mTimestepsPde0_1[model_solver_combo]).first;
                                    found = true;
                                    if ((mTimestepsPde0_1[model_solver_combo]).second)
                                    {
                                        accurate_enough = true;
                                    }
                                }
                            }
                            else if (pde_timestep==0.01)
                            {
                                std::map<std::pair<std::string, Solvers::Value>, std::pair<double,bool> >::iterator it = mTimestepsPde0_01.find(model_solver_combo);
                                if (it != mTimestepsPde0_01.end())
                                {
                                    ode_timestep = (mTimestepsPde0_01[model_solver_combo]).first;
                                    found = true;
                                    if ((mTimestepsPde0_01[model_solver_combo]).second)
                                    {
                                        accurate_enough = true;
                                    }
                                }
                            }
                            else
                            {
                                EXCEPTION("Summat went wrong. pde_timestep = " << pde_timestep << ", which wasn't in my files.");
                            }

                            if (!found)
                            {
                                WARNING("No suggested timestep found for " << model << " with '" << CellModelUtilities::GetSolverName(solver)
                                        << "' for pde timestep " << pde_timestep << ".\nSkipping this.");
                                continue;
                            }

                            if (!accurate_enough)
                            {
                                WARNING("Timestep " << ode_timestep << " for " << model << " with '" << CellModelUtilities::GetSolverName(solver) << "' for pde timestep " <<
                                        pde_timestep << " is not accurate enough, took longer than CVODE in calculate timesteps.\nSo skipping this.");
                                continue;
                            }
                        }
                        std::cout << "Model: " << model << " is being solved with " << CellModelUtilities::GetSolverName(solver)
                            << lookup_tables_description << " with ODE timestep = " << ode_timestep << "ms." << std::endl;
                        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(ode_timestep, pde_timestep, 0.1);

                        std::stringstream output_subfolder;
                        output_subfolder << "/results_pde_" <<  pde_timestep;

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

                        double total_elapsed_time;
                        double ode_elapsed_time;
                        try
                        {
                            Timer::Reset();
                            monodomain_problem.Solve();
                            total_elapsed_time = Timer::GetElapsedTime();
                            ode_elapsed_time = HeartEventHandler::GetElapsedTime(HeartEventHandler::SOLVE_ODES)/1000.0; // convert to seconds
                        }
                        catch (Exception& e)
                        {
                            std::cout << model << " failed to solve with ODE timestep " << ode_timestep << ", got: " << e.GetMessage() << std::endl << std::flush;
                            WARNING(model << " failed to solve with ODE timestep " << ode_timestep << ", got: " << e.GetMessage());
                            continue;
                        }
                        HeartEventHandler::Headings();
                        HeartEventHandler::Report();

                        if (PetscTools::AmMaster())
                        {
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

                                std::cout << "Model: " << model << " and PDE step " << pde_timestep << ":\nTime taken = " <<
                                        ode_elapsed_time << "/" << total_elapsed_time << " (ode/total)\nSquare error = " <<
                                        errors[0] << ", APD90 error = " << errors[1] << ", APD50 error = " << errors[2] <<
                                        ", APD30 error = " << errors[3] << ", V_max error = " << errors[4] << ", V_min error = " <<
                                        errors[5] << ", dVdt_max error = " << errors[6] << ", MRMS error = " << errors[7] << std::endl;

                                // Write to file too.
                                *p_file << model << "\t" << solver << "\t" << use_lookup_tables << "\t" << pde_timestep << "\t" << ode_timestep  << "\t"
                                        << total_elapsed_time << "\t" << ode_elapsed_time;
                                for (unsigned i=0; i<errors.size(); i++)
                                {
                                    *p_file << "\t" << errors[i];
                                }
                                *p_file << std::endl;
                            }
                            catch(Exception &e)
                            {
                                WARNING("Model " << model << ": analysis of voltage at last node failed.");
                            }
                        }
                    }

                    // Free memory for lookup tables if used
                    cell_factory.FreeLookupTableMemory();
                }
            }
        }

        if (PetscTools::AmMaster())
        {
            p_file->close();
        }
    }

private:

    /**
     * A map between the model/solver/pde-timestep and the ODE timestep/whether it's a 'refined enough' result.
     * See TestMonodomainCalculateRequiredTimestepsLiteratePaper.hpp
     *
     * For PDE step 0.1ms.
     */
    std::map<std::pair<std::string, Solvers::Value>, std::pair<double,bool> > mTimestepsPde0_1;

    /**
     * A map between the model/solver/pde-timestep and the ODE timestep/whether it's a 'refined enough' result.
     * See TestMonodomainCalculateRequiredTimestepsLiteratePaper.hpp
     *
     * For PDE step 0.01ms.
     */
    std::map<std::pair<std::string, Solvers::Value>, std::pair<double,bool> > mTimestepsPde0_01;

    /**
     * A helper method that populates mTimesteps from the stored data file in
     * Frontiers2014/test/data/required_steps.txt
     */
    void LoadTimestepFile()
    {
        FileFinder this_file(__FILE__);
        FileFinder summary_file("data/required_steps_tissue.txt", this_file);

        std::ifstream indata; // indata is like cin
        indata.open(summary_file.GetAbsolutePath().c_str()); // opens the file
        if(!indata.good())
        { // file couldn't be opened
            EXCEPTION("Couldn't open data file: " + summary_file.GetAbsolutePath());
        }

        while (indata.good())
        {
           std::string this_line;
           getline(indata, this_line);

           if (this_line=="" || this_line=="\r")
           {
               if (indata.eof())
               {    // If the blank line is the last line carry on OK.
                   break;
               }
               else
               {
                   EXCEPTION("No data found on this line");
               }
           }
           std::stringstream line(this_line);

           // Load a standard data line.
           std::string model_name;
           double tmp, pde_timestep, ode_timestep;
           bool is_satisfactory;
           int solver_index;

           line >> model_name;
           line >> pde_timestep;
           line >> solver_index;
           line >> ode_timestep;
           for (unsigned i=0; i<CellModelUtilities::NUM_ERROR_METRICS; i++)
           {
               line >> tmp;
           }
           line >> is_satisfactory;

           //std::cout << model_name << "\t" << solver_index << "\t" << pde_timestep << "\t" << ode_timestep << "\t" << is_satisfactory << "\n";

           Solvers::Value solver = (Solvers::Value)(solver_index); // We can read an int from this
           std::pair<std::string, Solvers::Value> model_solver_combo(model_name,solver);
           std::pair<double, bool> timestep_pair(ode_timestep, is_satisfactory);
           if (pde_timestep==0.1)
           {
               mTimestepsPde0_1[model_solver_combo] = timestep_pair;
           }
           else
           {
               mTimestepsPde0_01[model_solver_combo] = timestep_pair;
           }
        }

        if (!indata.eof())
        {
            EXCEPTION("A file reading error occurred");
        }

//        std::cout << "Other format:\n";
//        // Print to screen just to check they are correct...
//        for (std::map<std::pair<std::string, Solvers::Value>, std::pair<double, bool> >::iterator it = mTimestepsPde0_1.begin();
//             it!=mTimestepsPde0_1.end();
//             ++it)
//        {
//            if (!((it->second).second))
//            {
//                WARNING("Using non-satisfactory timestep for model " << (it->first).first << " and solver "
//                        << CellModelUtilities::GetSolverName((it->first).second));
//            }
//            // Print model, solver, timestep and whether it was satisfactory for a converged answer.
//            std::cout << (it->first).first << "\t'" << CellModelUtilities::GetSolverName((it->first).second) << "'\t" << (it->second).first << "\t" << (it->second).second << std::endl;
//        }

        std::cout << std::flush;
    }
};

#endif // TESTMONODOMAINSOLVINGTIMESLITERATEPAPER_HPP_
