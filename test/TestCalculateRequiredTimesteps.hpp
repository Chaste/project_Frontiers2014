/*

Copyright (c) 2005-2013, University of Oxford.
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

#include <cxxtest/TestSuite.h>

#include <fstream>
#include <algorithm>
#include <map>
#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#include <boost/pointer_cast.hpp>

#include "CellModelUtilities.hpp"

#include "FileFinder.hpp"
#include "AbstractCardiacCell.hpp"
#include "PetscTools.hpp"
#include "Warnings.hpp"
#include "IsNan.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestCalculateRequiredTimesteps : public CxxTest::TestSuite
{
public:
    /**
     * This test simulates each model with each (single-cell) method, reducing the time step used until the error metric
     * is within 5% of that achieved with CVODE in TestGeneratingReferenceData.hpp.
     */
    void TestCalculateTimesteps() throw (Exception)
    {
        /* First load up the error reference data from CVODE at different tolerances. */
        LoadErrorSummaryFile();

        /* Now iterate over model / solver combinations to find a suitable time step. */
        std::vector<FileFinder> models = CellModelUtilities::GetListOfModels();
        std::vector<Solvers::Value> solvers = boost::assign::list_of(Solvers::FORWARD_EULER)(Solvers::RUNGE_KUTTA_2)(Solvers::RUNGE_KUTTA_4)
                (Solvers::RUSH_LARSEN)(Solvers::GENERALISED_RUSH_LARSEN_1)(Solvers::GENERALISED_RUSH_LARSEN_2);

        // The result data structure
        std::map<std::pair<std::string, Solvers::Value>, std::pair<double,bool> > required_timesteps;

        // Create the output folder structure before isolating processes, to avoid race conditions
        OutputFileHandler test_base_handler("Frontiers/CalculateTimesteps/", false);
        BOOST_FOREACH(FileFinder& r_model, models)
        {
            OutputFileHandler model_handler("Frontiers/CalculateTimesteps/" + r_model.GetLeafNameNoExtension(), false);
        }

        PetscTools::IsolateProcesses();
        unsigned iteration = 0u;
        BOOST_FOREACH(FileFinder& r_model, models)
        {
            std::string model_name = r_model.GetLeafNameNoExtension();

            // Check if there is a reference trace for this model, and skip it if not
            if (mErrorResults.find(model_name) == mErrorResults.end())
            {
                continue;
            }

            BOOST_FOREACH(Solvers::Value solver, solvers)
            {
                if (iteration++ % PetscTools::GetNumProcs() != PetscTools::GetMyRank())
                {
                    continue; // Let another process do this combination
                }

                /* Generate the cell model from CellML. */
                std::stringstream folder_name;
                folder_name << "Frontiers/CalculateTimesteps/" << model_name << "/" << solver;
                OutputFileHandler handler(folder_name.str());
                boost::shared_ptr<AbstractCardiacCellInterface> p_cell = CellModelUtilities::CreateCellModel(r_model, handler, solver, false);
                double period = CellModelUtilities::GetDefaultPeriod(p_cell, 1000.0);

                /* Set up solver parameters. */
                double timestep = 0.2; // immediately divided by two so then a divisor of printing step.
                bool within_tolerance = false;
                std::vector<double> initial_conditions = p_cell->GetStdVecStateVariables();

                do
                {
                    timestep /= 2.0;
                    p_cell->SetTimestep(timestep);
                    p_cell->SetStateVariables(initial_conditions);

                    OdeSolution solution;
                    try
                    {
                        solution = p_cell->Compute(0.0, period, 0.1);
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

                    std::stringstream filename;
                    filename << model_name << "_deltaT_" << timestep;
                    solution.WriteToFile(handler.GetRelativePath(), filename.str(), "ms", 1, false, 16, false);
                    double error = CellModelUtilities::GetError(solution, model_name);
                    std::cout << model_name << ": '" << CellModelUtilities::GetSolverName(solver) << " error with timestep of " << timestep << "' = "
                            << error << ", required error = " << mErrorResults[model_name] << "." << std::endl;
                    within_tolerance = (error <= mErrorResults[model_name] * 1.05);
                }
                while (!within_tolerance && timestep >= 5e-6);

                // Record the timestep we needed
                std::pair<std::string, Solvers::Value> key(model_name, solver);
                std::pair<double, bool> result(timestep, within_tolerance);
                required_timesteps[key] = result;
            }
        }

        /* Write the summary file per process. */
        FileFinder this_file(__FILE__);
        FileFinder data_folder("data", this_file);
        out_stream p_file = test_base_handler.OpenOutputFile("required_steps_", PetscTools::GetMyRank(), ".txt");
        typedef std::pair<std::pair<std::string, Solvers::Value>, std::pair<double,bool> > ItemType;
        BOOST_FOREACH(ItemType item, required_timesteps)
        {
            // Writes model, solver, timestep, within_tolerance (tab separated)
            *p_file << item.first.first << "\t" << item.first.second << "\t" << item.second.first << "\t" << item.second.second << std::endl;
        }
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
                process_file_name << "required_steps_" << i << ".txt";
                std::ifstream process_file(process_file_name.str(), std::ios::binary);
                *p_combined_file << process_file.rdbuf();
            }
            p_combined_file->close();
            // Copy to repository
            FileFinder summary_file = test_base_handler.FindFile("required_steps.txt");
            summary_file.CopyTo(data_folder);
        }
    }

private:
    std::map<std::string, double> mErrorResults;

    /**
     * A helper method that populates mErrorResults from the stored data file in
     * Frontiers2014/data/reference_traces/error_summary.txt
     */
    void LoadErrorSummaryFile()
    {
        FileFinder this_file(__FILE__);
        FileFinder summary_file("data/reference_traces/error_summary.txt", this_file);

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
           double error_value;
           line >> model_name;
           line >> error_value;
           mErrorResults[model_name] = error_value;
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
