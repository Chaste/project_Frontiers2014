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

#ifndef TESTGENERATINGREFERENCEDATA_HPP_
#define TESTGENERATINGREFERENCEDATA_HPP_

/*
 * = Generate accurate reference traces and record a target error metric =
 *
 * On this wiki page we describe in detail the code that is used to run this example from the paper.
 *
 * For each cell model, this runs a single action potential at high fidelity (CVODE with tight tolerance and small max timestep)
 * as a 'gold standard' result for later comparison.
 *
 * Also, it runs a single action potential with looser tolerances as a "typical simulation" to give an accuracy benchmark used
 * to set time steps for other solvers in CalculateRequiredTimesteps.
 *
 * At the end of this test this information is copied into the file `test/data/error_summary.txt` within the project.
 *
 * == Code overview ==
 *
 * The first thing to do is to include the necessary header files.
 */

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

/* Next, we have the class that contains the code to run.
 * As with most papers using Chaste, this is written as a test suite within our test framework.
 * The `TestGenerateTraces` method does the work.
 */
class TestGeneratingReferenceData : public CxxTest::TestSuite
{
public:
    void TestGenerateTraces() throw (Exception)
    {
        /* Set up some references to various files and folders. */
        FileFinder this_file(__FILE__);
        // Where to copy reference traces within this project
        FileFinder repo_data("data/reference_traces", this_file);
        // Base folder for simulations to write output to
        OutputFileHandler test_base_handler("Frontiers/ReferenceTraces/", false);

        /* `CellModelUtilities` is a class of utility functions specifically for this project.
         * This one gets a list of all the CellML models included.
         */
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
            OutputFileHandler handler(test_base_handler.FindFile(model_name));
            boost::shared_ptr<AbstractCardiacCellInterface> p_cell = CellModelUtilities::CreateCellModel(r_model, handler, Solvers::CVODE_NUMERICAL_J, false);
            boost::shared_ptr<AbstractCvodeCell> p_cvode_cell = boost::dynamic_pointer_cast<AbstractCvodeCell>(p_cell);

            // If there's an analytic Jacobian available, we'll switch it on.
            if (p_cvode_cell->HasAnalyticJacobian())
            {
                p_cvode_cell->ForceUseOfNumericalJacobian(false);
            }

            /* Set up solver parameters. */
            p_cvode_cell->SetTolerances(1e-7 /* relative */, 1e-9 /* absolute */);

            /* Create a reference solution with high tolerances, and fine output (sampling every 0.1ms).
             * Note that the sampling interval is also the maximum time step CVODE is allowed to use.
             */
            double period = CellModelUtilities::GetDefaultPeriod(p_cell);
            OdeSolution solution;
            try
            {
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
                std::vector<std::pair<std::string, double> > properties;

                /* Calculate some summary statistics of the AP that was produced. */
                properties.push_back(std::pair<std::string, double>("Num_ODEs", (double)p_cvode_cell->GetNumberOfStateVariables()));
                properties.push_back(std::pair<std::string, double>("APD90", props.GetLastActionPotentialDuration(90.0)));
                properties.push_back(std::pair<std::string, double>("APD50", props.GetLastActionPotentialDuration(50.0)));
                properties.push_back(std::pair<std::string, double>("APD30", props.GetLastActionPotentialDuration(30.0)));
                properties.push_back(std::pair<std::string, double>("V_max", props.GetLastPeakPotential()));
                properties.push_back(std::pair<std::string, double>("V_min", props.GetLastRestingPotential()));
                properties.push_back(std::pair<std::string, double>("dVdt_max", props.GetLastMaxUpstrokeVelocity()));

                /* Save these to a dedicated file for this model, and copy to reference data folder in the repository. */
                out_stream p_summary_file = handler.OpenOutputFile(model_name + ".summary");
                for (unsigned i=0; i<properties.size(); i++)
                {
                    std::cout << properties[i].first  << " = " << properties[i].second << std::endl;
                    *p_summary_file << properties[i].first << "\t" << properties[i].second << std::endl;
                }
                p_summary_file->close();
                FileFinder summary_info = handler.FindFile(model_name + ".summary");
                summary_info.CopyTo(repo_data);
            }
            catch (const Exception& r_e)
            {
                WARNING("Action potential properties calculation failed for model " << model_name);
                continue;
            }
        }
    }
};

#endif // TESTGENERATINGREFERENCEDATA_HPP_
