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

#include "MathsCustomFunctions.hpp"
#include "ColumnDataReader.hpp"

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
    if (useLookupTables && p_cell->GetLookupTableCollection() == NULL)
    {
        EXCEPTION("No lookup tables implemented in optimised cell model " << rModelFile.GetLeafNameNoExtension());
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
    boost::shared_ptr<RegularStimulus> p_stim(new RegularStimulus(-25.5,5,1000,10));
    if (!p_cell->HasCellMLDefaultStimulus())
    {
        p_cell->SetStimulusFunction(p_stim);
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

double CellModelUtilities::GetError(const OdeSolution& rSolution, const std::string& rModelName)
{
    FileFinder this_file(__FILE__);
    FileFinder reference_folder("../test/data/reference_traces", this_file);

    ColumnDataReader data_reader(reference_folder, rModelName);
    std::vector<double> valid_times = data_reader.GetValues("Time");
    std::vector<double> valid_voltages = data_reader.GetValues("membrane_voltage");

    const std::vector<double>& r_new_times = rSolution.rGetTimes();
    std::vector<double> new_voltages = rSolution.GetAnyVariable("membrane_voltage");

    EXCEPT_IF_NOT(valid_times.size() == r_new_times.size());
    EXCEPT_IF_NOT(valid_voltages.size() == new_voltages.size());

    double error = 0.0;
    for (unsigned i=0; i<valid_times.size(); i++)
    {
        EXCEPT_IF_NOT(CompareDoubles::WithinAbsoluteTolerance(valid_times[i], r_new_times[i], 1e-12));
        error += SmallPow((valid_voltages[i] - new_voltages[i]), 2u);
    }
    return error;
}
