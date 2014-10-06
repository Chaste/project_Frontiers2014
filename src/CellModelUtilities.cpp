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
    FileFinder this_file(__FILE__);
    FileFinder model_folder("../cellml", this_file);
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

    // This method needs to work with a direct pointer to a cell, so that methods in DynamicModelCellFactory can use it.
    CellModelUtilities::SetCellModelSolver(p_cell.get(), solver);

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


void CellModelUtilities::SetCellModelSolver(AbstractCardiacCellInterface* pCell,
                                            Solvers::Value solver)
{
    // Check that we have the features we're expecting, and specify solver if not built-in
    switch (solver)
    {
        case Solvers::CVODE_ANALYTIC_J:
        case Solvers::CVODE_NUMERICAL_J:
        {
            AbstractCvodeCell* p_cvode_cell = dynamic_cast<AbstractCvodeCell*>(pCell);
            assert(p_cvode_cell != NULL);
            if (solver == Solvers::CVODE_ANALYTIC_J)
            {
                if (!p_cvode_cell->GetUseAnalyticJacobian())
                {
                    EXCEPTION("No analytic Jacobian available for cell model " << p_cvode_cell->GetSystemName());
                }
            }
            else
            {
                p_cvode_cell->ForceUseOfNumericalJacobian();
            }
            p_cvode_cell->SetTolerances(1e-4 /* relative */, 1e-6 /* absolute */);
            p_cvode_cell->SetMaxSteps(10000000);
            boost::shared_ptr<RegularStimulus> p_reg_stim = boost::dynamic_pointer_cast<RegularStimulus>(pCell->GetStimulusFunction());
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
            pCell->SetSolver(p_rk2_solver);
            break;
        }
        case Solvers::RUNGE_KUTTA_4:
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_rk4_solver(new RungeKutta4IvpOdeSolver());
            pCell->SetSolver(p_rk4_solver);
            break;
        }
        default:
            break; // Nothing to do here
    }
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
            WARNING("Setting start time for stimulus to 0ms as was " << p_reg_stim->GetStartTime() << "ms which is greater than period " << p_reg_stim->GetPeriod()
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

void CellModelUtilities::SetCvodeTolerances(AbstractCardiacCellInterface* pCell, unsigned index)
{
    if (index >= 7u)
    {
        /*
         * Since 1e-7, 1e-9 was used to generate the reference traces, we should never need to go below this,
         * so throw an Exception if this happens.
         */
        EXCEPTION("SetCvodeTolerances: Asking for more refinement of CVODE than we think should be necessary!");
    }

    std::vector<std::pair<double,double> > possible_tolerances;
    for (double relative_tol = 1e-2; relative_tol > 1e-8; relative_tol/=10.0)
    {
        std::pair<double, double> tolerance_pair(relative_tol, relative_tol/100.0);
        possible_tolerances.push_back(tolerance_pair);
    }

    AbstractCvodeCell* p_cvode_cell = dynamic_cast<AbstractCvodeCell*>(pCell);

    if (p_cvode_cell)
    {
        p_cvode_cell->SetTolerances(possible_tolerances[index].first, possible_tolerances[index].second);
    }
    else
    {
        EXCEPTION("You need to pass a pointer to (something that can be cast to) an AbstractCvodeCell to the method SetCvodeTolerances().");
    }

}

std::vector<double> CellModelUtilities::GetErrors(const OdeSolution& rSolution, const std::string& rModelName)
{
    FileFinder this_file(__FILE__);
    FileFinder reference_folder("../test/data/reference_traces", this_file);

    ColumnDataReader data_reader(reference_folder, rModelName);
    std::vector<double> valid_times = data_reader.GetValues("Time");
    std::vector<double> valid_voltages = data_reader.GetValues("membrane_voltage");

    const std::vector<double>& r_new_times = rSolution.rGetTimes();
    std::vector<double> new_voltages = rSolution.GetAnyVariable("membrane_voltage");

    return GetCommonErrorCalculations(valid_times,valid_voltages,r_new_times,new_voltages);
}

std::vector<double> CellModelUtilities::GetTissueErrors(const std::vector<double>& rTestTimes,
                                                        const std::vector<double>& rTestVoltages,
                                                        const std::string& rModelName,
                                                        const double& rPdeTimestep)
{

    FileFinder this_file(__FILE__);
    FileFinder reference_folder("../test/data/reference_traces", this_file);

    // Load up the CVODE with resetting as the definitive trace.
    std::stringstream filename_stream;
    filename_stream << rModelName << "_tissue_pde_" << rPdeTimestep << "_h_0.01_Reset.dat";

    // Load the reference traces from file
    FileFinder reference_trace(filename_stream.str(), reference_folder);
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

    return GetCommonErrorCalculations(valid_times,valid_voltages,rTestTimes,rTestVoltages);
}

std::vector<double> CellModelUtilities::GetCommonErrorCalculations(const std::vector<double>& rValidTimes,
                                                                   const std::vector<double>& rValidVoltages,
                                                                   const std::vector<double>& rTestTimes,
                                                                   const std::vector<double>& rTestVoltages)
{
    std::vector<double> errors;

    EXCEPT_IF_NOT(rValidTimes.size() == rTestTimes.size());
    EXCEPT_IF_NOT(rValidVoltages.size() == rTestVoltages.size());
    EXCEPT_IF_NOT(rValidVoltages.size() == rValidTimes.size());

    double square_error = 0.0;
    double mixed_root_mean_square = 0.0;
    for (unsigned i=0; i<rValidTimes.size(); i++)
    {
        EXCEPT_IF_NOT(CompareDoubles::WithinAbsoluteTolerance(rValidTimes[i], rTestTimes[i], 1e-12));
        double tmp = SmallPow((rValidVoltages[i] - rTestVoltages[i]), 2u);
        square_error += tmp;
        mixed_root_mean_square += tmp / SmallPow(1 + fabs(rValidVoltages[i]), 2u);
    }
    mixed_root_mean_square = sqrt(mixed_root_mean_square/(double)(rValidTimes.size()));

    errors.push_back(square_error);

    CellProperties reference_properties(rValidVoltages, rValidTimes);
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

    errors.push_back(mixed_root_mean_square);

    // N.B. If you add any more error metrics, be sure to increase NUM_ERROR_METRICS in the .hpp file.

    return errors;
}
