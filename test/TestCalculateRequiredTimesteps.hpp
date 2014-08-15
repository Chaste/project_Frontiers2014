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
#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#include <boost/pointer_cast.hpp>

#include "FileFinder.hpp"
#include "CellModelUtilities.hpp"
#include "AbstractCardiacCell.hpp"
#include "Warnings.hpp"

class TestCalculateRequiredTimesteps : public CxxTest::TestSuite
{
public:
    void TestCalculateTimesteps() throw (Exception)
    {
        /**
         * First load up the error reference data from CVODE at different tolerances
         */
        LoadErrorSummaryFile();

        std::vector<FileFinder> models = CellModelUtilities::GetListOfModels();

        // Loop over each model we have
        BOOST_FOREACH(FileFinder& r_model, models)
        {
            /* Generate the cell model from CellML. */
            std::string model_name = r_model.GetLeafNameNoExtension();

            // Check if there is a reference trace for this model...
            bool reference_traces_exist = false;
            for (std::map<std::string,double>::iterator it = mErrorResults.begin();
                 it!=mErrorResults.end();
                 ++it)
            {
                if (it->first == model_name)
                {
                    reference_traces_exist = true;
                    break;
                }
            }
            // If there isn't, then skip it.
            if (!reference_traces_exist)
            {
                continue;
            }

            OutputFileHandler handler("Frontiers/CalculateTimesteps/" + model_name);
            std::vector<std::string> options = boost::assign::list_of("--Wu"); // Just a 'plain' model.
            boost::shared_ptr<AbstractCardiacCellInterface> p_cell = CellModelUtilities::CreateCellModel(r_model, handler, options);
            double period = CellModelUtilities::GetDefaultPeriod(p_cell, 1000.0);

            /* Set up solver parameters. */

            double timestep = 0.2; // immediately divided by two so then a divisor of printing step.
            bool within_tolerance = false;
            std::vector<double> initial_conditions = p_cell->GetStdVecStateVariables();

            do
            {
                timestep/=2.0;
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
                std::cout << model_name << ": 'Forward Euler error with timestep of " << timestep << "' = "
                        << error << ", required error = " << mErrorResults[model_name] << "." << std::endl << std::flush;
                within_tolerance = (error<=mErrorResults[model_name]);
            }
            while (!within_tolerance && timestep>=5e-6);

            // Record the timestep we needed
            std::pair<double, bool> result(timestep, within_tolerance);
            mRequiredForwardEulerTimestep[model_name] = result;
        }

        /* Write the summary file, and copy it to the repository */
        FileFinder this_file(__FILE__);
        FileFinder data_folder("data", this_file);
        OutputFileHandler handler("Frontiers/CalculateTimesteps", false);
        out_stream p_file = handler.OpenOutputFile("forward_euler_steps.txt");
        for (std::map<std::string,std::pair<double,bool> >::iterator it = mRequiredForwardEulerTimestep.begin();
             it!=mRequiredForwardEulerTimestep.end();
             ++it)
        {
            *p_file << it->first << "\t" << (it->second).first << "\t" << (it->second).second << std::endl;
        }
        p_file->close();

        FileFinder summary_file = handler.FindFile("forward_euler_steps.txt");
        summary_file.CopyTo(data_folder);
    }

private:
    std::map<std::string, double> mErrorResults;
    std::map<std::string, std::pair<double,bool> > mRequiredForwardEulerTimestep; // timestep and whether it met CVODE tolerances.

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
