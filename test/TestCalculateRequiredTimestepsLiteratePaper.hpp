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

#ifndef TESTCALCULATEREQUIREDTIMESTEPS_HPP_
#define TESTCALCULATEREQUIREDTIMESTEPS_HPP_

/*
 * = Calculate the required ODE solver timesteps to meet target accuracy =
 *
 * On this wiki page we describe in detail the code that is used to run this example from the paper.
 *
 * This test simulates each model with each (single-cell) numerical method,
 * reducing the time step used until the error metric is within 5% of that achieved
 * with CVODE using slack tolerances in PaperTutorials/Frontiers2014/GeneratingReferenceData
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

// This header is needed to allow us to run in parallel
#include "PetscSetupAndFinalize.hpp"

class TestCalculateRequiredTimesteps : public CxxTest::TestSuite
{
public:
    void TestCalculateTimesteps() throw (Exception)
    {
        /* First load up the error reference data from CVODE at different tolerances. */
        LoadErrorSummaryFile();

        /* The model / solver combinations to find a suitable time step for. */
        std::vector<FileFinder> models = CellModelUtilities::GetListOfModels();
        std::vector<Solvers::Value> solvers = boost::assign::list_of(Solvers::FORWARD_EULER)(Solvers::RUNGE_KUTTA_2)
                                                                    (Solvers::RUNGE_KUTTA_4)(Solvers::RUSH_LARSEN)
                                                                    (Solvers::GENERALISED_RUSH_LARSEN_1)
                                                                    (Solvers::GENERALISED_RUSH_LARSEN_2);

        /* Create the output folder structure before isolating processes, to avoid race conditions. */
        OutputFileHandler test_base_handler("Frontiers/CalculateTimesteps/", false);
        BOOST_FOREACH(FileFinder& r_model, models)
        {
            OutputFileHandler model_handler("Frontiers/CalculateTimesteps/" + r_model.GetLeafNameNoExtension(), false);
        }

        /* Each process writes its results to a separate file. */
        FileFinder this_file(__FILE__);
        FileFinder data_folder("data", this_file);
        out_stream p_file = test_base_handler.OpenOutputFile("required_steps_", PetscTools::GetMyRank(), ".txt");
        /* Ensure time steps are written at high precision */
        *p_file << std::setiosflags(std::ios::scientific) << std::setprecision(16);

        /* Iterate over model/solver combinations, distributed over processes. */
        PetscTools::IsolateProcesses();
        unsigned iteration = 0u;
        BOOST_FOREACH(FileFinder& r_model, models)
        {
            std::string model_name = r_model.GetLeafNameNoExtension();

            // Check if there is a reference trace for this model, and skip it if not
            if (mErrorResults.find(model_name) == mErrorResults.end())
            {
                std::cout << "No reference error metric for model " << model_name << "; skipping." << std::endl;
                WARNING("No reference error metric for model " << model_name << "; skipping.");
                continue;
            }

            /* Iterate over each available solver, using a handy boost method */
            BOOST_FOREACH(Solvers::Value solver, solvers)
            {
                /* This simple if allows us to parallelise the sweep to speed up running it */
                if (iteration++ % PetscTools::GetNumProcs() != PetscTools::GetMyRank())
                {
                    continue; // Let another process do this combination
                }
                std::cout << "Calculating timestep for " << model_name << " with solver " << CellModelUtilities::GetSolverName(solver) << std::endl;

                /* Generate the cell model from CellML. */
                std::stringstream folder_name;
                folder_name << "Frontiers/CalculateTimesteps/" << model_name << "/" << solver;
                OutputFileHandler handler(folder_name.str());
                boost::shared_ptr<AbstractCardiacCellInterface> p_cell = CellModelUtilities::CreateCellModel(r_model, handler, solver, false);
                double period = CellModelUtilities::GetDefaultPeriod(p_cell);

                /* Set up solver parameters. */
                double sampling_time = 0.1;
                unsigned timestep_divisor = 1u;
                bool within_tolerance = false;
                std::vector<double> initial_conditions = p_cell->GetStdVecStateVariables();

                /* Keep solving with smaller timesteps until we meet the desired accuracy */
                do
                {
                    double timestep = sampling_time / timestep_divisor;
                    timestep_divisor *= 2; // Ready for next iteration
                    p_cell->SetTimestep(timestep);
                    /* Make sure the model is in exactly the same state at the beginning of every run */
                    p_cell->SetStateVariables(initial_conditions);

                    OdeSolution solution;
                    /* We enclose the solve in a try-catch, as large timesteps can lead to solver crashes */
                    try
                    {
                        solution = p_cell->Compute(0.0, period, sampling_time);
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
                        std::cout << model_name << " failed to solve with timestep " << timestep << ", got: " << e.GetMessage() << std::endl << std::flush;
                        continue;
                    }

                    /* Calculate the error associated with this simulated trace (compare to reference) */
                    std::stringstream filename;
                    filename << model_name << "_deltaT_" << timestep;
                    solution.WriteToFile(handler.GetRelativePath(), filename.str(), "ms", 1, false, 16, false);
                    try
                    {
                        std::vector<double> errors = CellModelUtilities::GetError(solution, model_name);
                        std::cout << model_name << ": '" << CellModelUtilities::GetSolverName(solver) << " error with timestep of " << timestep << "' = "
                                << errors[0] << ", required error = " << mErrorResults[model_name][0] << "." << std::endl;
                        within_tolerance = (errors[0] <= mErrorResults[model_name][0] * 1.05);

                        /* Write result to file, tab separated */
                        *p_file << model_name << "\t" << solver << "\t" << timestep;
                        for (unsigned i=0; i<errors.size(); i++)
                        {
                            *p_file << "\t" << errors[i];
                        }
                        *p_file << "\t" << within_tolerance << std::endl;
                    }
                    catch (const Exception& r_e)
                    {
                        std::cout << "Exception computing error for model " << model_name << " solver " << CellModelUtilities::GetSolverName(solver);
                        std::cout << r_e.GetMessage() << std::endl;
                        WARNING(r_e.GetMessage());
                    }
                }
                while (!within_tolerance && timestep_divisor <= 2048);
                /* The above means we allow at most 12 refinements of the time step (2^12^ = 2048). */
            }
        }

        /* Close each process' results file. */
        p_file->close();

        /* Turn off process isolation and wait for all files to be written. */
        PetscTools::IsolateProcesses(false);
        PetscTools::Barrier("TestCalculateTimesteps");

        /* Master process writes the concatenated file. */
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
            /* Copy to repository for storage and use by other tests */
            FileFinder summary_file = test_base_handler.FindFile("required_steps.txt");
            summary_file.CopyTo(data_folder);
        }
    }

private:
    /**
     * A map between the model name and the square error associated with a 'loose' CVODE solve,
     * as generated by TestGeneratingReferenceData.hpp
     */
    std::map<std::string, std::vector<double> > mErrorResults;

    /**
     * A helper method that populates mErrorResults from the stored data file in
     * Frontiers2014/test/data/error_summary.txt
     */
    void LoadErrorSummaryFile()
    {
        FileFinder this_file(__FILE__);
        FileFinder summary_file("data/error_summary.txt", this_file);

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
           std::vector<double> error_values(7);
           line >> model_name;
           for (unsigned i=0; i<7 ; i++)
           {
               line >> error_values[i];
           }
           mErrorResults[model_name] = error_values;
        }

        if (!indata.eof())
        {
            EXCEPTION("A file reading error occurred");
        }

// Print to screen just to check they are correct...
//        for (std::map<std::string,double>::iterator it = mErrorResults.begin();
//             it!=mErrorResults.end();
//             ++it)
//        {
//            std::cout << it->first << "\t" << it->second << std::endl;
//        }
    }
};

#endif // TESTCALCULATEREQUIREDTIMESTEPS_HPP_
