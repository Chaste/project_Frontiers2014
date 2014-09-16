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

#ifndef CELLMODELUTILITIES_HPP_
#define CELLMODELUTILITIES_HPP_

#include <vector>
#include <string>

#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractCardiacCellInterface.hpp"

/**
 * Structure encapsulating the enumeration of cell model solvers.
 * This allows us to write things like Solvers::CVODE for readability.
 */
struct Solvers
{
    /**
     * What solvers are available for cardiac cell models.
     */
    enum Value
    {
        CVODE_ANALYTIC_J,          /**< CVODE with analytically calculated Jacobian */
        CVODE_NUMERICAL_J,         /**< CVODE with numerical approximation to Jacobian */
        FORWARD_EULER,             /**< Forward Euler */
        BACKWARD_EULER,            /**< Tissue simulation specific Backward Euler */
        RUNGE_KUTTA_2,             /**< 2nd order Runge-Kutta */
        RUNGE_KUTTA_4,             /**< 4th order Runge-Kutta */
        RUSH_LARSEN,               /**< Rush-Larsen method */
        GENERALISED_RUSH_LARSEN_1, /**< Generalised Rush-Larsen method */
        GENERALISED_RUSH_LARSEN_2  /**< Generalised Rush-Larsen method */
    };
};

/**
 * This class defines various utility methods used within this project, for tasks such as loading CellML models
 * and generating code with particular solver/optimisation properties.
 */
class CellModelUtilities
{
private:
    /**
     * Dynamically convert a CellML file into C++ code, compile it and load the resulting model object.
     *
     * @param rModelFile  the CellML file to convert
     * @param rOutputDir  handler for the folder in which to create generated code
     * @param rOptions  extra options for the code generation tool
     *
     * @return a cardiac cell for Chaste to use.
     */
    static boost::shared_ptr<AbstractCardiacCellInterface> CreateCellModel(const FileFinder& rModelFile,
                                                                           OutputFileHandler& rOutputDir,
                                                                           const std::vector<std::string>& rOptions);
public:
    /**
     * Read the CellML project folder and return all the .cellml files defined therein.
     */
    static std::vector<FileFinder> GetListOfModels();

    /**
     * Get the human readable name associated with a solver code.
     * @param solver  the solver code
     * @return  the solver's name
     */
    static std::string GetSolverName(Solvers::Value solver);

    /**
     * Dynamically convert a CellML file into C++ code, compile it and load the resulting model object.
     *
     * @param rModelFile  the CellML file to convert
     * @param rOutputDir  handler for the folder in which to create generated code
     * @param solver  which cell model solver to use
     * @param useLookupTables  whether to use lookup tables to speed up simulations
     *
     * @return a cardiac cell for Chaste to use.
     */
    static boost::shared_ptr<AbstractCardiacCellInterface> CreateCellModel(const FileFinder& rModelFile,
                                                                           OutputFileHandler& rOutputDir,
                                                                           Solvers::Value solver,
                                                                           bool useLookupTables);

    /**
     * Get the pacing cycle length to use when simulating a cell.
     * If the cell has no RegularStimulus attached, we use a default of 1000ms.
     * If the stimulus period is not a multiple of 1ms, we round down to the nearest millisecond,
     * to ensure that time steps or sampling intervals used for simulation will divide the simulation duration exactly.
     *
     * @param pCell  the cell to extract stimulus information from
     * @return The default period
     */
    static double GetDefaultPeriod(boost::shared_ptr<AbstractCardiacCellInterface> pCell);

    /**
     * Compare a simulation result against reference data, and compute error metrics.
     * We use:
     *  * Sum of square error in each voltage time point,
     *  * absolute error in:
     *    > APD90,
     *    > APD50,
     *    > APD30,
     *    > V_max,
     *    > V_min, (perhaps always initial condition, so not that useful)
     *    > dVdt_max.
     *
     * Note that a lot of these may end up being zero, due to output time steps making the calculations
     * quite coarse.
     *
     * @param rSolution  the new simulation results
     * @param rModelName  the name of the model simulated, used to find the reference results in the project folder
     *
     * @return the error metrics listed above.
     */
    static std::vector<double> GetErrors(const OdeSolution& rSolution,
                                         const std::string& rModelName);

    /**
     * Compare a simulation result against reference data (the AP at the last node),
     * and compute error metrics.
     *
     * We use:
     *  * Sum of square error in each voltage time point,
     *  * absolute error in:
     *    > APD90,
     *    > APD50,
     *    > APD30,
     *    > V_max,
     *    > V_min, (perhaps always initial condition, so not that useful)
     *    > dVdt_max.
     *
     * @param rTestTimes  the test simulation printing times
     * @param rTestVoltages  the test simulation voltages at these times on last node.
     * @param rModelName  the name of the model simulated, used to find the reference results in the project folder
     * @param rPdeTimestep  the PDE timestep we are using, so we load the same reference trace to avoid confusing PDE and ODE convergence.
     *
     * @return the error metrics listed above.
     */
    static std::vector<double> GetTissueErrors(const std::vector<double>& rTestTimes,
                                               const std::vector<double>& rTestVoltages,
                                               const std::string& rModelName,
                                               const double& rPdeTimestep);

    /**
     * Return suitable error metrics (the ones listed above in GetError documentation), for
     * this model for a slack vs. tight CVODE solve.
     *
     * @param tissue  Whether we want to load the monodomain reference traces (defaults to false)
     * @param requestedPdeTimestep  The PDE timestep to load if tissue==true.
     *
     * @return a map between the model name, and the vector of error metrics
     */
    static std::map<std::string,  std::vector<double> > LoadErrorSummaryFile(bool tissue=false,
                                                                             double requestedPdeTimestep=0.01);

    /**
     * Set the CVODE tolerances from a table that looks like
     *
     * Index    RelTol  AbsTol
     * 1        1e-2    1e-4
     * 2        1e-3    1e-5
     * 3        1e-4    1e-6
     * 4        1e-5    1e-7
     * ...      ...     ...
     *
     * @param pCell  A pointer to something that can be cast as an AbstractCvodeCell (will throw if not)
     * @param index  The index of the row of the table to use (i.e. how many steps of refinement to use).
     */
    static void SetCvodeTolerances(boost::shared_ptr<AbstractCardiacCellInterface> pCell,
                                   unsigned index);
};

#endif // CELLMODELUTILITIES_HPP_
