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

#include "CellModelUtilities.hpp"

#include <algorithm>
#include <cmath>
#include <boost/foreach.hpp>

#include "MathsCustomFunctions.hpp"
#include "ColumnDataReader.hpp"
#include "Warnings.hpp"

#include "CellProperties.hpp"
#include "CellMLLoader.hpp"
#include "AbstractUntemplatedParameterisedSystem.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

bool finder_sort(const FileFinder& rf1, const FileFinder& rf2)
{
    return rf1.GetAbsolutePath() < rf2.GetAbsolutePath();
}

std::vector<FileFinder> CellModelUtilities::GetListOfModels()
{
    FileFinder model_folder("projects/Frontiers2014/cellml", RelativeTo::ChasteSourceRoot);
    std::vector<FileFinder> models = model_folder.FindMatches("*.cellml");
    std::sort(models.begin(), models.end(), finder_sort);
    return models;
}

std::string CellModelUtilities::GetSolverName(Solvers::Value solver)
{
    std::string name;
    switch (solver)
    {
        case Solvers::CVODE_ANALYTIC_J:
            name = "CVODE (analytic Jacobian)";
            break;
        case Solvers::CVODE_NUMERICAL_J:
            name = "CVODE (numerical Jacobian)";
            break;
        case Solvers::BACKWARD_EULER:
            name = "Backward Euler";
            break;
        case Solvers::RUSH_LARSEN:
            name = "Rush-Larsen";
            break;
        case Solvers::GENERALISED_RUSH_LARSEN_1:
            name = "Generalised Rush-Larsen 1";
            break;
        case Solvers::GENERALISED_RUSH_LARSEN_2:
            name = "Generalised Rush-Larsen 2";
            break;
        case Solvers::FORWARD_EULER:
            name = "Forward Euler";
            break;
        case Solvers::RUNGE_KUTTA_2:
            name = "Runge-Kutta (2nd order)";
            break;
        case Solvers::RUNGE_KUTTA_4:
            name = "Runge-Kutta (4th order)";
            break;
        default:
            EXCEPTION("Unknown solver code!");
            break;
    }
    return name;
}

boost::shared_ptr<AbstractCardiacCellInterface> CellModelUtilities::CreateCellModel(const FileFinder& rModelFile,
                                                                                    OutputFileHandler& rOutputDir,
                                                                                    Solvers::Value solver,
                                                                                    bool useLookupTables)
{
    // Figure out code generation options
    std::vector<std::string> options;
    switch (solver)
    {
        case Solvers::CVODE_ANALYTIC_J:
        case Solvers::CVODE_NUMERICAL_J:
            options.push_back("--cvode");
            break;
        case Solvers::BACKWARD_EULER:
            options.push_back("--backward-euler");
            break;
        case Solvers::RUSH_LARSEN:
            options.push_back("--rush-larsen");
            break;
        case Solvers::GENERALISED_RUSH_LARSEN_1:
            options.push_back("--grl1");
            break;
        case Solvers::GENERALISED_RUSH_LARSEN_2:
            options.push_back("--grl2");
            break;
        case Solvers::FORWARD_EULER:
        case Solvers::RUNGE_KUTTA_2:
        case Solvers::RUNGE_KUTTA_4:
            break;
        default:
            EXCEPTION("Unexpected solver passed to CreateCellModel.");
            break;
    }
    if (useLookupTables)
    {
        options.push_back("--opt");
    }

    // Load the cell model
    boost::shared_ptr<AbstractCardiacCellInterface> p_cell = CreateCellModel(rModelFile, rOutputDir, options);

    // Check that we have the features we're expecting, and specify solver if not built-in
    switch (solver)
    {
        case Solvers::CVODE_ANALYTIC_J:
        case Solvers::CVODE_NUMERICAL_J:
        {
            boost::shared_ptr<AbstractCvodeCell> p_cvode_cell = boost::dynamic_pointer_cast<AbstractCvodeCell>(p_cell);
            if (solver == Solvers::CVODE_ANALYTIC_J)
            {
                if (!p_cvode_cell->GetUseAnalyticJacobian())
                {
                    EXCEPTION("No analytic Jacobian available for cell model " << rModelFile.GetLeafNameNoExtension());
                }
            }
            else
            {
                p_cvode_cell->ForceUseOfNumericalJacobian();
            }
            p_cvode_cell->SetTolerances(1e-4 /* relative */, 1e-6 /* absolute */);
            p_cvode_cell->SetMaxSteps(10000000);
            boost::shared_ptr<RegularStimulus> p_reg_stim = boost::dynamic_pointer_cast<RegularStimulus>(p_cell->GetStimulusFunction());
            // If the max time step for cvode is not set then it defaults to the sampling timestep,
            // if this is large then CVODE can 'skip' the stimulus and never notice it, so should generally
            // be set to the stimulus duration.
            if (p_reg_stim)
            {
                p_cvode_cell->SetMaxTimestep(p_reg_stim->GetDuration());
            }
            break;
        }
        case Solvers::RUNGE_KUTTA_2:
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_rk2_solver(new RungeKutta2IvpOdeSolver());
            p_cell->SetSolver(p_rk2_solver);
            break;
        }
        case Solvers::RUNGE_KUTTA_4:
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_rk4_solver(new RungeKutta4IvpOdeSolver());
            p_cell->SetSolver(p_rk4_solver);
            break;
        }
        default:
            break; // Nothing to do here
    }
    // Check that the lookup tables exist if they should
    if (useLookupTables)
    {
        AbstractLookupTableCollection* p_tables = p_cell->GetLookupTableCollection();
        if (p_tables == NULL || p_tables->GetKeyingVariableNames().empty())
        {
            EXCEPTION("No lookup tables implemented in optimised cell model " << rModelFile.GetLeafNameNoExtension());
        }
        BOOST_FOREACH(std::string keying_var, p_tables->GetKeyingVariableNames())
        {
            if (p_tables->GetNumberOfTables(keying_var) == 0u)
            {
                // This should never happen (the code generation shouldn't allow it) but just in case...
                EXCEPTION("No lookup tables for '" << keying_var << "' implemented in optimised cell model " << rModelFile.GetLeafNameNoExtension());
            }
        }
    }

    return p_cell;
}


boost::shared_ptr<AbstractCardiacCellInterface> CellModelUtilities::CreateCellModel(
        const FileFinder& rModelFile,
        OutputFileHandler& rOutputDir,
        const std::vector<std::string>& rOptions)
{
    // Add --Wu to options if missing, and test if we're using CVODE or not
    std::vector<std::string> args(rOptions);
    if (std::find(args.begin(), args.end(), "--Wu") == args.end())
    {
        args.push_back("--Wu");
    }
    bool use_cvode = (std::find(args.begin(), args.end(), "--cvode") != args.end());

    // Load the model, with its default stimulus
    CellMLLoader loader(rModelFile, rOutputDir, args);
    boost::shared_ptr<AbstractCardiacCellInterface> p_cell;
    if (use_cvode)
    {
        p_cell = loader.LoadCvodeCell();
    }
    else
    {
        p_cell = loader.LoadCardiacCell();
    }

    // Make a default stimulus to work with
    // (if the model has a stimulus current, but no defaults tagged, it was ending up with an empty stimulus (e.g. HH1952)).
    // TODO Put in a check for this before using a Cell model where the RHS includes the stimulus current.
    // TODO Consider if generated code should expose the pycml:is-self-excitatory annotation, or an indication that dV/dt doesn't use the stimulus.
    if (!p_cell->HasCellMLDefaultStimulus())
    {
        boost::shared_ptr<RegularStimulus> p_stim(new RegularStimulus(-25.5,5,1000,10));
        p_cell->SetStimulusFunction(p_stim);
    }
    // Check for unusual stimulus patterns where the first stimulus occurs after a full period
    boost::shared_ptr<RegularStimulus> p_reg_stim = boost::dynamic_pointer_cast<RegularStimulus>(p_cell->GetStimulusFunction());
    if (p_reg_stim)
    {
        if (p_reg_stim->GetStartTime() > p_reg_stim->GetPeriod())
        {
            WARNING("Setting start time for stimulus to 0 as was greater than period " << p_reg_stim->GetPeriod()
                    << " in model " << rModelFile.GetLeafNameNoExtension() << ".");
            p_reg_stim->SetStartTime(0.0);
        }
    }

    return p_cell;
}

double CellModelUtilities::GetDefaultPeriod(boost::shared_ptr<AbstractCardiacCellInterface> pCell)
{
    double period = 1000.0;
    boost::shared_ptr<RegularStimulus> p_stim = boost::dynamic_pointer_cast<RegularStimulus>(pCell->GetStimulusFunction());
    if (p_stim)
    {
        period = p_stim->GetPeriod();
    }
    else
    {
        std::string name = boost::dynamic_pointer_cast<AbstractUntemplatedParameterisedSystem>(pCell)->GetSystemName();
        std::cout << "Note: no regular stimulus found for " << name << " so using default period of 1000ms." << std::endl;
    }
    if (floor(period) != period)
    {
        std::string name = boost::dynamic_pointer_cast<AbstractUntemplatedParameterisedSystem>(pCell)->GetSystemName();
        std::cout << "Note: rounding period for " << name << " from " << period << " to " << floor(period) << "." << std::endl;
        period = floor(period);
    }
    return period;
}

std::vector<double> CellModelUtilities::GetErrors(const OdeSolution& rSolution, const std::string& rModelName)
{
    std::vector<double> errors;
    FileFinder this_file(__FILE__);
    FileFinder reference_folder("../test/data/reference_traces", this_file);

    ColumnDataReader data_reader(reference_folder, rModelName);
    std::vector<double> valid_times = data_reader.GetValues("Time");
    std::vector<double> valid_voltages = data_reader.GetValues("membrane_voltage");

    const std::vector<double>& r_new_times = rSolution.rGetTimes();
    std::vector<double> new_voltages = rSolution.GetAnyVariable("membrane_voltage");

    EXCEPT_IF_NOT(valid_times.size() == r_new_times.size());
    EXCEPT_IF_NOT(valid_voltages.size() == new_voltages.size());

    double square_error = 0.0;
    for (unsigned i=0; i<valid_times.size(); i++)
    {
        EXCEPT_IF_NOT(CompareDoubles::WithinAbsoluteTolerance(valid_times[i], r_new_times[i], 1e-12));
        square_error += SmallPow((valid_voltages[i] - new_voltages[i]), 2u);
    }
    errors.push_back(square_error);

    CellProperties reference_properties(valid_voltages, valid_times);
    CellProperties test_properties(new_voltages, r_new_times);

    errors.push_back(reference_properties.GetLastActionPotentialDuration(90.0)
                     - test_properties.GetLastActionPotentialDuration(90.0));

    errors.push_back(reference_properties.GetLastActionPotentialDuration(50.0)
                     - test_properties.GetLastActionPotentialDuration(50.0));

    errors.push_back(reference_properties.GetLastActionPotentialDuration(30.0)
                     - test_properties.GetLastActionPotentialDuration(30.0));

    errors.push_back(reference_properties.GetLastPeakPotential()
                     - test_properties.GetLastPeakPotential());

    errors.push_back(reference_properties.GetLastRestingPotential()
                     - test_properties.GetLastRestingPotential());

    errors.push_back(reference_properties.GetLastMaxUpstrokeVelocity()
                     - test_properties.GetLastMaxUpstrokeVelocity());

    return errors;
}

std::vector<double> CellModelUtilities::GetTissueErrors(const std::vector<double>& rTestTimes,
                                                        const std::vector<double>& rTestVoltages,
                                                        const std::string& rModelName)
{
    std::vector<double> errors;
    FileFinder this_file(__FILE__);
    FileFinder reference_folder("../test/data/reference_traces", this_file);

    // Load the reference traces from file
    FileFinder reference_trace(rModelName + "_tissue.dat", reference_folder);
    EXCEPT_IF_NOT(reference_trace.IsFile());

    std::vector<double> valid_times;
    std::vector<double> valid_voltages;
    {
        std::ifstream indata; // indata is like cin
        indata.open(reference_trace.GetAbsolutePath().c_str()); // opens the file
        if(!indata.good())
        { // file couldn't be opened
            EXCEPTION("Couldn't open data file: " + reference_trace.GetAbsolutePath());
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
           double time;
           double voltage;
           line >> time;
           line >> voltage;
           valid_times.push_back(time);
           valid_voltages.push_back(voltage);
        }

        if (!indata.eof())
        {
            EXCEPTION("A file reading error occurred");
        }

        EXCEPT_IF_NOT(valid_times.size() == rTestTimes.size());
        EXCEPT_IF_NOT(valid_voltages.size() == rTestVoltages.size());
    }

    double square_error = 0.0;
    for (unsigned i=0; i<valid_times.size(); i++)
    {
        EXCEPT_IF_NOT(CompareDoubles::WithinAbsoluteTolerance(valid_times[i], rTestTimes[i], 1e-12));
        square_error += SmallPow((valid_voltages[i] - rTestVoltages[i]), 2u);
    }
    errors.push_back(square_error);

    CellProperties reference_properties(valid_voltages, valid_times);
    CellProperties test_properties(rTestVoltages, rTestTimes);

    errors.push_back(reference_properties.GetLastActionPotentialDuration(90.0)
                     - test_properties.GetLastActionPotentialDuration(90.0));

    errors.push_back(reference_properties.GetLastActionPotentialDuration(50.0)
                     - test_properties.GetLastActionPotentialDuration(50.0));

    errors.push_back(reference_properties.GetLastActionPotentialDuration(30.0)
                     - test_properties.GetLastActionPotentialDuration(30.0));

    errors.push_back(reference_properties.GetLastPeakPotential()
                     - test_properties.GetLastPeakPotential());

    errors.push_back(reference_properties.GetLastRestingPotential()
                     - test_properties.GetLastRestingPotential());

    errors.push_back(reference_properties.GetLastMaxUpstrokeVelocity()
                     - test_properties.GetLastMaxUpstrokeVelocity());

    return errors;
}

/**
 * A helper method that populates an error results structure from the stored data file in
 * Frontiers2014/test/data/error_summary.txt
 */
std::map<std::string, std::vector<double> > CellModelUtilities::LoadErrorSummaryFile()
{
    std::map<std::string, std::vector<double> > error_results;

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
       error_results[model_name] = error_values;
    }

    if (!indata.eof())
    {
        EXCEPTION("A file reading error occurred");
    }

    return error_results;
}


