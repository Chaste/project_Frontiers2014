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

#ifndef TESTGENERATINGREFERENCEDATA_HPP_
#define TESTGENERATINGREFERENCEDATA_HPP_

// The testing framework we use
#include <cxxtest/TestSuite.h>

// Standard C++ libraries
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>

// Headers specific to this project
#include "CellModelUtilities.hpp"

// General Chaste headers
#include "AbstractCvodeCell.hpp"
#include "CellProperties.hpp"
#include "PetscTools.hpp"
#include "Warnings.hpp"

// This header is needed to allow us to run in parallel
#include "PetscSetupAndFinalize.hpp"

/**
 * For each cell model, runs a single action potential at high fidelity (CVODE with tight tolerance and small max timestep)
 * as a 'gold standard' result for later comparison.
 *
 * Also, runs a single action potential with looser tolerances as a "typical simulation" to give an accuracy benchmark used
 * to set time steps for other solvers in TestCalculateRequiredTimesteps.hpp. At the end of this test this information is
 * copied into test/data/reference_traces/error_summary.txt
 */
class TestGeneratingReferenceData : public CxxTest::TestSuite
{
public:
    void TestGenerateTraces() throw (Exception)
    {
        FileFinder this_file(__FILE__);
        FileFinder repo_data("data/reference_traces", this_file);
        OutputFileHandler test_base_handler("Frontiers/ReferenceTraces/", false);
        std::map<std::string, double> error_results;

        std::vector<FileFinder> models = CellModelUtilities::GetListOfModels();

        /* Iterate over the available models, handling each one on a separate process if running in parallel. */
        PetscTools::IsolateProcesses();

        for (unsigned i=0; i<models.size(); ++i)
        {
            if (i % PetscTools::GetNumProcs() != PetscTools::GetMyRank())
            {
                continue; // Let another process handle this model
            }

            /* Generate the cell model from CellML. */
            FileFinder& r_model = models[i];
            std::string model_name = r_model.GetLeafNameNoExtension();
            OutputFileHandler handler("Frontiers/ReferenceTraces/" + model_name);
            std::vector<std::string> options = boost::assign::list_of("--Wu")("--cvode");
            boost::shared_ptr<AbstractCardiacCellInterface> p_cell = CellModelUtilities::CreateCellModel(r_model, handler, options);
            double period = CellModelUtilities::GetDefaultPeriod(p_cell);

            /* Set up solver parameters. */
            boost::shared_ptr<AbstractCvodeCell> p_cvode_cell = boost::dynamic_pointer_cast<AbstractCvodeCell>(p_cell);
            p_cvode_cell->SetTolerances(1e-7 /* relative */, 1e-9 /* absolute */);

            /* Create a reference solution with high tolerances, and fine output (sampling every 0.1ms). */
            OdeSolution solution;
            try
            {
                // Note that the sampling interval of 0.1ms is also the maximum time step CVODE is allowed to use
                solution = p_cvode_cell->Compute(0.0, period, 0.1);
            }
            catch (const Exception &e)
            {
                WARNING(model_name << " model failed to solve. It gave the error:\n" << e.GetMessage());
                continue;
            }

            /* Write solution to file, and copy to repository folder. */
            solution.WriteToFile(handler.GetRelativePath(), model_name, "ms", 1, false, 16, false);
            FileFinder results_dat = handler.FindFile(model_name + ".dat");
            results_dat.CopyTo(repo_data);
            FileFinder results_info = handler.FindFile(model_name + ".info");
            results_info.CopyTo(repo_data);

            /* Check that the solution looks like an action potential. */
            try
            {
                std::vector<double> voltages = solution.GetVariableAtIndex(p_cell->GetVoltageIndex());
                CellProperties props(voltages, solution.rGetTimes());
                props.GetLastActionPotentialDuration(90.0);
            }
            catch (const Exception& r_e)
            {
                WARNING("No action potential produced for model " << model_name);
                continue;
            }

            /* Perform another simulation with lower CVODE tolerances, to give us an error bound to set other solvers' time steps with. */
            p_cvode_cell->SetTolerances(1e-4 /* relative */, 1e-6 /* absolute */);
            p_cvode_cell->ResetToInitialConditions();
            solution = p_cvode_cell->Compute(0.0, period, 0.1);

            /* Write this solution to file so we can compare graphs, and print the error metric. */
            solution.WriteToFile(handler.GetRelativePath(), model_name + "_rough", "ms", 1, false, 16, false);
            double error = CellModelUtilities::GetError(solution, model_name);
            std::cout << "Model " << model_name << " square error " << error << std::endl;
            error_results[model_name] = error;
        }

        /* Write a file containing the target error metric for each model.
         * Since we may have run simulations in parallel, each process writes its own file, and the master processes then
         * concatenates these and copies the resulting single file into the Chaste repository.
         */

        // Each process writes its own file, at high precision
        out_stream p_error_file = test_base_handler.OpenOutputFile("error_summary_", PetscTools::GetMyRank(), ".txt");
        *p_error_file << std::setiosflags(std::ios::scientific) << std::setprecision(16);
        typedef std::pair<std::string, double> StringDoublePair;
        BOOST_FOREACH(StringDoublePair error, error_results)
        {
            *p_error_file << error.first << "\t" << error.second << std::endl;
        }
        p_error_file->close();

        // Turn off process isolation and wait for all files to be written
        PetscTools::IsolateProcesses(false);
        PetscTools::Barrier("TestGenerateTraces");

        // Master process writes the concatenated file
        if (PetscTools::AmMaster())
        {
            out_stream p_combined_file = test_base_handler.OpenOutputFile("error_summary.txt", std::ios::out | std::ios::trunc | std::ios::binary);
            for (unsigned i=0; i<PetscTools::GetNumProcs(); ++i)
            {
                std::stringstream process_file_name;
                process_file_name << test_base_handler.GetOutputDirectoryFullPath() << "error_summary_" << i << ".txt";
                std::ifstream process_file(process_file_name.str().c_str(), std::ios::binary);
                TS_ASSERT(process_file.is_open());
                TS_ASSERT(process_file.good());
                *p_combined_file << process_file.rdbuf();
            }
            p_combined_file->close();
            // Copy to repository
            FileFinder error_summary_file = test_base_handler.FindFile("error_summary.txt");
            error_summary_file.CopyTo(repo_data);
        }
    }
};

#endif // TESTGENERATINGREFERENCEDATA_HPP_
