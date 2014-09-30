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

#ifndef TESTODESOLVINGTIMES_HPP_
#define TESTODESOLVINGTIMES_HPP_

/*
 * = Benchmark ODE solving times with different numerical methods =
 *
 * On this wiki page we describe in detail the code that is used to run this example from the paper.
 *
 * For each cell model, we load the required timestep for the required target accuracy,
 * as calculated by PaperTutorials/Frontiers2014/CalculateRequiredTimesteps
 *
 * We then time solving all models under most solvers, for 10 paces.
 *
 * == Code overview ==
 *
 * The first thing to do is to include the necessary header files.
 */

// The testing framework we use
#include <cxxtest/TestSuite.h>

// Standard C++ libraries
#include <iostream>
#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

// Headers specific to this project
#include "CellModelUtilities.hpp"

// General Chaste headers
#include "Timer.hpp"
#include "Warnings.hpp"
#include "ChasteBuildRoot.hpp"
#include "AbstractCvodeCell.hpp"

// This header is needed to allow us to run in parallel
#include "PetscSetupAndFinalize.hpp"

static unsigned NUM_RUNS = 3u;

class TestOdeSolvingTimes : public CxxTest::TestSuite
{
    /* These are some type aliases to save typing. */
    typedef boost::tuple<std::string, Solvers::Value, bool> KeyType; // Keys in the map types below
    typedef std::map<KeyType, double> TimestepMapType; // Type of the map giving suggested timesteps
    typedef std::map<KeyType, bool> TimestepOkMapType; // Type of the map indicating convergence solutions

public:
    /* This is a small visual check that we have the expected CellML models available. */
	void TestListingModels() throw (Exception)
	{
		std::vector<FileFinder> models(CellModelUtilities::GetListOfModels());
		std::cout << "Available models:" << std::endl;
		BOOST_FOREACH(FileFinder& r_model, models)
		{
		    std::cout << "  " << r_model.GetLeafNameNoExtension() << std::endl;
		}
		TS_ASSERT_LESS_THAN(20u, models.size());
	}

	/* This is the main test that actually does the benchmarking. */
	void TestSolvingTimes() throw (Exception)
	{
	    /* First load the suggested time steps to use. */
        LoadTimestepFile();

        const double required_mrms_error = 0.05; // 5%

        /* The model / solver combinations to find a suitable time step for. */
        std::vector<FileFinder> models = CellModelUtilities::GetListOfModels();
        std::vector<Solvers::Value> solvers = boost::assign::list_of
                (Solvers::CVODE_ANALYTIC_J)(Solvers::CVODE_NUMERICAL_J)
                (Solvers::FORWARD_EULER)(Solvers::BACKWARD_EULER)
                (Solvers::RUNGE_KUTTA_2)(Solvers::RUNGE_KUTTA_4)(Solvers::RUSH_LARSEN)
                (Solvers::GENERALISED_RUSH_LARSEN_1)(Solvers::GENERALISED_RUSH_LARSEN_2);

        /* Create the output folder structure before isolating processes, to avoid race conditions. */
        OutputFileHandler test_base_handler("Frontiers/SingleCellTimings/", false);
        BOOST_FOREACH(FileFinder& r_model, models)
        {
            OutputFileHandler model_handler("Frontiers/SingleCellTimings/" + r_model.GetLeafNameNoExtension(), false);
        }

        /* Each process writes its timings to a separate file as we go along, in case of catastrophe. */
        out_stream p_file = test_base_handler.OpenOutputFile(ChasteBuildType() + "_timings_", PetscTools::GetMyRank(), ".txt");
        *p_file << std::setiosflags(std::ios::scientific) << std::setprecision(8);

        /* Each process writes its catastrophies too. */
        out_stream p_errors_file = test_base_handler.OpenOutputFile(ChasteBuildType() + "_errors_", PetscTools::GetMyRank(), ".txt");
        *p_errors_file << std::setiosflags(std::ios::scientific) << std::setprecision(8);

        /* Iterate over model/solver combinations, distributed over processes. */
        PetscTools::IsolateProcesses();
        unsigned iteration = 0u;
        BOOST_FOREACH(FileFinder& r_model, models)
        {
            std::string model_name = r_model.GetLeafNameNoExtension();

            BOOST_FOREACH(Solvers::Value solver, solvers)
            {
                if (iteration++ % PetscTools::GetNumProcs() != PetscTools::GetMyRank())
                {
                    continue; // Let another process do this combination
                }

                std::string solver_name = CellModelUtilities::GetSolverName(solver);

                bool cvode_solver = ((solver==Solvers::CVODE_ANALYTIC_J) || (solver==Solvers::CVODE_NUMERICAL_J));

                std::vector<bool> lookup_table_options = boost::assign::list_of(false)(true);
                BOOST_FOREACH(bool use_lookup_tables, lookup_table_options)
                {
                    std::string using_tables = (use_lookup_tables ? " and lookup tables" : "");

                    /* Get timestep to use if available.
                     * Note that the CreateCellModel method below sets suitable tolerances for CVODE.
                     */
                    double suggested_timestep = 0.0; // Signifies 'not set'
                    double time_to_simulate = 10; // seconds, a default
                    KeyType key(model_name, solver, use_lookup_tables);
                    TimestepMapType::iterator dt_iter = mTimesteps.find(key);
                    if (dt_iter != mTimesteps.end())
                    {
                        suggested_timestep = dt_iter->second;
                        time_to_simulate = mSimulationTime[key];
                        std::cout << "Simulating " << model_name << " with solver '" << solver_name << "'" << using_tables
                                  << " using timestep " << suggested_timestep << std::endl;
                    }
                    else
                    {
                        std::cout << "No timestep result found for " << model_name << " with solver '" << solver_name << "'" << using_tables
                                  << ", skipping it." << std::endl;
                        continue;
                    }

                    try
                    {
                        /* Generate the cell model from CellML. */
                        std::stringstream folder_name;
                        folder_name << "Frontiers/SingleCellTimings/" << model_name << "/" << solver;
                        if (use_lookup_tables)
                        {
                            folder_name << "_Opt";
                        }
                        OutputFileHandler handler(folder_name.str());
                        boost::shared_ptr<AbstractCardiacCellInterface> p_cell = CellModelUtilities::CreateCellModel(r_model, handler, solver, use_lookup_tables);
                        double period = CellModelUtilities::GetDefaultPeriod(p_cell);

                        if (cvode_solver)
                        {
                            // We hijacked the timestep in the log file to provide a refinement level.
                            unsigned refinement_idx = (unsigned)(suggested_timestep);
                            CellModelUtilities::SetCvodeTolerances(p_cell.get(), refinement_idx);
                        }
                        else
                        {
                            p_cell->SetTimestep(suggested_timestep);
                        }

                        /*
                         * Run a single pace to check accuracy.
                         * This will also set up lookup tables if they are requested, and prevent the
                         * one-off setup cost being timed
                         */
                        OdeSolution solution = p_cell->Compute(0.0, period, 0.1);
                        solution.WriteToFile(handler.GetRelativePath(), model_name, "ms", 1, false, 16, false);

                        /* Just double check the accuracy is what we expect before a timing run */
                        std::vector<double> errors = CellModelUtilities::GetErrors(solution, model_name);
                        std::cout << "Model " << model_name << " solver '" << solver_name << "'" << using_tables << " MRMS error = " << errors[7] << std::endl;
                        if (errors[7] > required_mrms_error)
                        {
                            WARNING("Model " << model_name << " with solver '" << solver_name << "'" << using_tables
                                    << " did not reach error target. Wanted " << required_mrms_error << ", got " << errors[7] << ".");
                        }

                        /* Time simulating multiple paces. */
                        double elapsed_time = TimeSimulation(p_cell, time_to_simulate);
                        std::cout << "Model " << model_name << " solver '" << solver_name << "'" << using_tables << " took time " << elapsed_time << "s per simulated sec" << std::endl;

                        /* Record the result. */
                        *p_file << model_name << "\t" << solver << "\t" << use_lookup_tables << "\t" << elapsed_time;
                        for (unsigned i=0; i<errors.size(); i++)
                        {
                            *p_file << "\t" << errors[i];
                        }
                        *p_file << std::endl;

                        // Free memory for lookup tables, if used
                        if (use_lookup_tables)
                        {
                            p_cell->GetLookupTableCollection()->FreeMemory();
                        }
                    }
                    catch (const Exception& r_e)
                    {
                        std::stringstream error_message;
                        error_message << "Error simulating model " << model_name << " with solver '" << solver_name << "'" <<  using_tables << ": " << r_e.GetMessage();
                        WARNING(error_message.str());
                        std::cout << error_message.str() << std::endl;
                        *p_errors_file << error_message.str() << std::endl;
                    }
                }
            }
        }

        /* Close each process' results file. */
        *p_file << "# Complete" << std::endl;
        p_file->close();
        p_errors_file->close();

        /* Turn off process isolation and wait for all files to be written. */
        PetscTools::IsolateProcesses(false);
        PetscTools::Barrier("TestSolvingTimes");

        /* Master process writes the concatenated files. */
        if (PetscTools::AmMaster())
        {
            // Do both the timings and errors files.
            std::vector<std::string> file_types = boost::assign::list_of("timings")("errors");
            BOOST_FOREACH(std::string file_type, file_types)
            {
                out_stream p_combined_file = test_base_handler.OpenOutputFile(ChasteBuildType() + "_" + file_type + ".txt", std::ios::out | std::ios::trunc | std::ios::binary);
                for (unsigned i=0; i<PetscTools::GetNumProcs(); ++i)
                {
                    std::stringstream process_file_name;
                    process_file_name << test_base_handler.GetOutputDirectoryFullPath() << ChasteBuildType() << "_" << file_type << "_" << i << ".txt";
                    std::ifstream process_file(process_file_name.str().c_str(), std::ios::binary);
                    TS_ASSERT(process_file.is_open());
                    TS_ASSERT(process_file.good());
                    *p_combined_file << process_file.rdbuf();
                }
                p_combined_file->close();
            }
        }
	}

private:
	/*
	 * Utility methods used by the tests above go here.
	 */

	/**
	 * Find out how long it takes to simulate the given model.
	 * The cell will be reset to initial conditions prior to simulation.
	 * We assume the cell has a regular square wave stimulus defined.
	 *
	 * Note that we don't check the results are sensible, or do a pre-simulation to avoid
	 * counting lookup tables setup.
	 *
	 * We convert the result so that it tells us how many wall clock seconds are required
	 * for a simulation of one second of electrophysiology.
	 *
	 * @param pCell  the cell model to simulate, with solver attached
	 * @param timeToSimulate  the amount of time to simulate in seconds
	 * @return  number of wall clock seconds to simulate one second's worth of this model activity.
	 */
	double TimeSimulation(boost::shared_ptr<AbstractCardiacCellInterface> pCell,
	                      double timeToSimulate)
	{
	    // Drop any fractional part of the simulation time (converted to ms) so it is divided by the timestep
	    double millisecs_to_simulate = floor(timeToSimulate * 1000.0);

	    //std::cout << "Simulating " << millisecs_to_simulate << "ms of activity." << std::endl;

	    // Run NUM_RUNS simulations and record the quickest
	    double minimum = DBL_MAX;
	    for (unsigned i=0; i<NUM_RUNS; i++)
	    {
            boost::dynamic_pointer_cast<AbstractUntemplatedParameterisedSystem>(pCell)->ResetToInitialConditions();
            Timer::Reset();
            pCell->SolveAndUpdateState(0.0, millisecs_to_simulate);
            double elapsed_time = Timer::GetElapsedTime();
            if (elapsed_time < minimum)
            {
                minimum = elapsed_time;
            }
	    }

	    /* Convert the elapsed time into a time per simulated second. */
	    return minimum/timeToSimulate;
	}


    /**
     * A map between the model/solver/use of lookup tables, and the timestep which gave a 'refined enough' result.
     * See TestCalculateRequiredTimesteps.hpp
     */
    TimestepMapType mTimesteps;

    /**
     * A map between the model/solver and the number of paces we should simulate.
     * This is set by looking at how long one pace took in #LoadTimestepFile.
     */
    TimestepMapType mSimulationTime;

    /**
     * A map between the model/solver/use of look uptables, and whether the timestep in #mTimesteps
     * gave a 'refined enough' result.
     * See TestCalculateRequiredTimesteps.hpp
     */
    TimestepOkMapType mTimestepIsSatisfactory;

    /**
     * A helper method that populates mTimesteps from the stored data file in
     * Frontiers2014/test/data/required_steps.txt
     */
    void LoadTimestepFile()
    {
        FileFinder this_file(__FILE__);
        FileFinder summary_file("data/required_steps.txt", this_file);

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
           double tmp, timestep;
           bool is_satisfactory, lt_used;
           int solver_index;
           double time_taken;

           line >> model_name;
           line >> solver_index;
           line >> lt_used;
           line >> timestep;
           for (unsigned i=0; i<8; i++)
           {
               line >> tmp;
           }
           line >> is_satisfactory;
           line >> time_taken;

           Solvers::Value solver = (Solvers::Value)(solver_index); // We can read an int from this.
           KeyType key(model_name, solver, lt_used);
           mTimesteps[key] = timestep;
           mTimestepIsSatisfactory[key] = is_satisfactory; // We always have the most refined last in the file, so this overwriting should work.

           // Now calculate a sensible amount of time to 'time' each model/solver for.
           // We'll aim for 5 seconds of simulation
           double target_time = 5.0; // seconds
           mSimulationTime[key] = target_time/time_taken;
        }

        if (!indata.eof())
        {
            EXCEPTION("A file reading error occurred");
        }

        // Warn for any model/solver combination not satisfactory
        for (TimestepOkMapType::iterator it = mTimestepIsSatisfactory.begin();
             it!=mTimestepIsSatisfactory.end();
             ++it)
        {
            if (!it->second)
            {
                WARNING("Using non-satisfactory timestep for model " << (it->first).get<0>() << " and solver "
                        << CellModelUtilities::GetSolverName((it->first).get<1>())
                        << ((it->first).get<2>() ? " with lookup tables" : ""));
            }
        }

        // Print to screen just to check they are correct...
        for (TimestepMapType::iterator it = mTimesteps.begin();
             it!=mTimesteps.end();
             ++it)
        {
            // Print model, solver, timestep and whether it was satisfactory for a converged answer.
            std::cout << (it->first).get<0>() << "\t'" << CellModelUtilities::GetSolverName((it->first).get<1>()) << "'\t" << (it->first).get<2>() << "\t"
                      << it->second << "\t" << mTimestepIsSatisfactory[it->first] << std::endl;
        }
    }

};

#endif // TESTODESOLVINGTIMES_HPP_
