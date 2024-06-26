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

#ifndef TESTCALCULATEREQUIREDTIMESTEPS_HPP_
#define TESTCALCULATEREQUIREDTIMESTEPS_HPP_

/*
 * = Calculate the required ODE solver timesteps to meet target accuracy =
 *
 * On this wiki page we describe in detail the code that is used to run this example from the paper.
 *
 * This test simulates each model with each (single-cell) numerical method,
 * reducing the time step used until the error metric is within 5% of that achieved
 * with CVODE using slack tolerances in GeneratingReferenceData.
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

// General Chaste headers
#include "FileFinder.hpp"
#include "AbstractCardiacCell.hpp"
#include "PetscTools.hpp"
#include "Warnings.hpp"
#include "IsNan.hpp"
#include "Timer.hpp"

// This header is needed to allow us to run in parallel
#include "PetscSetupAndFinalize.hpp"

class TestCalculateRequiredTimesteps : public CxxTest::TestSuite
{
public:
    void TestCalculateTimesteps() throw (Exception)
    {
        const double required_mrms_error = 0.05; // 5%

        /* The model / solver combinations to find a suitable time step for. */
        std::vector<FileFinder> models = CellModelUtilities::GetListOfModels();
        std::vector<Solvers::Value> solvers = boost::assign::list_of
                (Solvers::CVODE_ANALYTIC_J)(Solvers::CVODE_NUMERICAL_J)
                (Solvers::FORWARD_EULER)(Solvers::BACKWARD_EULER)
                (Solvers::RUNGE_KUTTA_2)(Solvers::RUNGE_KUTTA_4)
                (Solvers::RUSH_LARSEN)(Solvers::GENERALISED_RUSH_LARSEN_1)(Solvers::GENERALISED_RUSH_LARSEN_2);

        /* Create the output folder structure before isolating processes, to avoid race conditions. */
        OutputFileHandler test_base_handler("Frontiers/CalculateTimesteps/", false);
        BOOST_FOREACH(FileFinder& r_model, models)
        {
            OutputFileHandler model_handler(test_base_handler.FindFile(r_model.GetLeafNameNoExtension()), false);
        }

        /* Each process writes its results to a separate file. */
        out_stream p_file = test_base_handler.OpenOutputFile("required_steps_", PetscTools::GetMyRank(), ".txt");
        // Ensure time steps are written at high precision
        *p_file << std::setiosflags(std::ios::scientific) << std::setprecision(16);

        /* Iterate over model/solver combinations, distributed over processes. */
        PetscTools::IsolateProcesses();
        unsigned iteration = 0u;
        BOOST_FOREACH(FileFinder& r_model, models)
        {
            std::string model_name = r_model.GetLeafNameNoExtension();

            /* We calculate separate timesteps for with and without lookup tables,
             * since for some model/solver combinations including lookup tables can push things outside the stable regime.
             */
            std::vector<bool> lookup_table_options = boost::assign::list_of(false)(true);
            BOOST_FOREACH(bool use_lookup_tables, lookup_table_options)
            {
                std::string using_tables = (use_lookup_tables ? " and lookup tables" : "");

                /* Iterate over each available solver.
                 * This is the inner loop since it iterates over 9 items, making it unlikely to be divisible by the number
                 * of processes and hence hopefully giving a more even distribution of work.
                 */
                BOOST_FOREACH(Solvers::Value solver, solvers)
                {
                    bool cvode_solver = ((solver==Solvers::CVODE_ANALYTIC_J) || (solver==Solvers::CVODE_NUMERICAL_J));

                    /* This simple if allows us to parallelise the sweep to speed up running it. */
                    if (iteration++ % PetscTools::GetNumProcs() != PetscTools::GetMyRank())
                    {
                        continue; // Let another process do this combination
                    }
                    std::cout << "Calculating timestep for " << model_name << " with solver " << CellModelUtilities::GetSolverName(solver) << using_tables << std::endl;

                    /* Generate the cell model from CellML. */
                    std::stringstream folder_name;
                    folder_name << model_name << "/" << solver;
                    if (use_lookup_tables)
                    {
                        folder_name << "_Opt";
                    }
                    OutputFileHandler handler(test_base_handler.FindFile(folder_name.str()));
                    boost::shared_ptr<AbstractCardiacCellInterface> p_cell;
                    try
                    {
                        p_cell = CellModelUtilities::CreateCellModel(r_model, handler, solver, use_lookup_tables);
                    }
                    catch (Exception &e)
                    {
                        // This probably happens if there is no Analytic Jacobian available.
                        // Move on to the next model/solver combination.
                        continue;
                    }

                    /* Set up solver parameters. */
                    double sampling_time = 0.1;
                    unsigned timestep_divisor = 1u;
                    unsigned refinement_idx = 1u;
                    bool within_tolerance = false;
                    std::vector<double> initial_conditions = p_cell->GetStdVecStateVariables();
                    double period = CellModelUtilities::GetDefaultPeriod(p_cell);

                    /* For the one step models: keep solving with smaller timesteps until we meet the desired accuracy.
                     * For CVODE, use the timestep divisor to index a vector of tolerance pairs. */
                    do
                    {
                        double timestep = sampling_time / timestep_divisor;
                        std::stringstream refinement_description_stream;
                        if (cvode_solver)
                        {
                            if (refinement_idx>=7u)
                            {
                                // CVODE is struggling, probably because it is hitting a singularity.
                                break;
                            }
                            CellModelUtilities::SetCvodeTolerances(p_cell.get(), refinement_idx);
                            timestep = (double)(refinement_idx); // We just hijack this variable for the log filenames.
                            refinement_description_stream << " CVODE tolerances of " << pow(10.0,-1.0-refinement_idx) <<
                                    ", " << pow(10.0,-3.0-refinement_idx);
                        }
                        else
                        {
                            p_cell->SetTimestep(timestep);
                            refinement_description_stream << " timestep of " << timestep << "ms";
                        }
                        std::string refinement_description = refinement_description_stream.str();

                        std::cout << "Computing error for model " << model_name << " solver '" << CellModelUtilities::GetSolverName(solver) << "'"
                                  << using_tables << " with" << refinement_description << std::endl;

                        /* Prepare the loop variables for the next iteration here, since a subsequent block contains a
                         * `continue` statement, which would skip the adjustment if it happened later.
                         */
                        timestep_divisor *= 2;
                        refinement_idx++;

                        /* Make sure the model is in exactly the same state at the beginning of every run. */
                        p_cell->SetStateVariables(initial_conditions);

                        OdeSolution solution;
                        double time_taken;
                        /* We enclose the solve in a try-catch, as large timesteps can lead to solver crashes. */
                        try
                        {
                            Timer::Reset();
                            solution = p_cell->Compute(0.0, period, sampling_time);
                            time_taken = Timer::GetElapsedTime();
                            // Check that the solution doesn't contain any non-numerical values, indicative of solver failure.
                            std::vector<double> voltages = solution.GetAnyVariable("membrane_voltage");
                            std::vector<double>& r_times = solution.rGetTimes();
                            for (unsigned i=0; i<voltages.size(); i++)
                            {
                                if (std::isnan(voltages[i]))
                                {
                                    EXCEPTION("Voltage at t(" << r_times[i] << ") = " << voltages[i] << ".");
                                }
                            }
                        }
                        catch (const Exception &e)
                        {
                            std::cout << model_name << " failed to solve with" << refinement_description << ".\n Got: " << e.GetMessage() << std::endl << std::flush;
                            continue;
                        }

                        /* Calculate the error associated with this simulated trace (compared to the reference trace). */
                        std::stringstream filename;
                        filename << model_name << "_deltaT_" << timestep;
                        solution.WriteToFile(handler.GetRelativePath(), filename.str(), "ms", 1, false, 16, false);
                        try
                        {
                            std::vector<double> errors = CellModelUtilities::GetErrors(solution, model_name);
                            assert(errors.size()==8u);
                            std::cout << "MRMS error = " << errors[7] << " (required = " << required_mrms_error << ")." << std::endl;
                            within_tolerance = (errors[7] <= required_mrms_error);

                            // Write result to file, tab separated
                            *p_file << model_name << "\t" << solver << "\t" << use_lookup_tables << "\t" << timestep;
                            for (unsigned i=0; i<errors.size(); i++)
                            {
                                *p_file << "\t" << errors[i];
                            }
                            *p_file << "\t" << within_tolerance << "\t" << time_taken << std::endl;
                        }
                        catch (const Exception& r_e)
                        {
                            std::cout << "Exception computing error for model " << model_name << " solver " << CellModelUtilities::GetSolverName(solver) << using_tables;
                            std::cout << r_e.GetMessage() << std::endl;
                            WARNING(r_e.GetMessage());
                        }
                    }
                    while (!within_tolerance && timestep_divisor <= 2048);
                    /* The above means we allow at most 12 refinements of the time step (2^12^ = 2048). */

                    /* Free memory for lookup tables, if used.
                     * To avoid excessive memory use in tissue simulations, the lookup tables for a particular cell model are
                     * stored in a separate singleton class, which can thus be shared by all cells being simulated in a
                     * particular process.  Since here we are instead simulating many different models in turn, this is not so
                     * ideal, as we end up with a lookup table collection for each one persisting until the program exits.  This
                     * could lead to memory exhaustion if a fine table granularity is used, so we tell Chaste to delete the
                     * tables when we know we've finished with them.
                     */
                    if (use_lookup_tables)
                    {
                        p_cell->GetLookupTableCollection()->FreeMemory();
                    }
                }
            }
        }

        /* Close each process' results file. */
        p_file->close();

        /* Turn off process isolation and wait for all files to be written. */
        PetscTools::IsolateProcesses(false);
        PetscTools::Barrier("TestCalculateTimesteps");

        /* Master process writes the concatenated file, and copies it into the project's data folder. */
        if (PetscTools::AmMaster())
        {
            out_stream p_combined_file = test_base_handler.OpenOutputFile("required_steps.txt", std::ios::out | std::ios::trunc | std::ios::binary);
            for (unsigned i=0; i<PetscTools::GetNumProcs(); ++i)
            {
                std::stringstream process_file_name;
                process_file_name << test_base_handler.GetOutputDirectoryFullPath() << "required_steps_" << i << ".txt";
                std::ifstream process_file(process_file_name.str().c_str(), std::ios::binary);
                TS_ASSERT(process_file.is_open());
                TS_ASSERT(process_file.good());
                *p_combined_file << process_file.rdbuf();
            }
            p_combined_file->close();
            /* Copy to repository for storage and use by other tests. */
            FileFinder this_file(__FILE__);
            FileFinder data_folder("data", this_file);
            FileFinder summary_file = test_base_handler.FindFile("required_steps.txt");
            summary_file.CopyTo(data_folder);
        }
    }
};

#endif // TESTCALCULATEREQUIREDTIMESTEPS_HPP_
