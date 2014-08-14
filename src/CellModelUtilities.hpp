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
 * This class defines various utility methods used within this project, for tasks such as loading CellML models
 * and generating code with particular solver/optimisation properties.
 */
class CellModelUtilities
{
public:
    /**
     * Read the CellML project folder and return all the .cellml files defined therein.
     */
    static std::vector<FileFinder> GetListOfModels();

    /**
     * Dynamically convert a CellML file into C++ code, compile it and load the resulting model object.
     *
     * @param rModelFile  the CellML file to convert
     * @param rOutputDir  handler for the folder in which to create generated code
     * @param rOptions  extra options for the code generation tool
     */
    static boost::shared_ptr<AbstractCardiacCellInterface> CreateCellModel(const FileFinder& rModelFile,
                                                                           OutputFileHandler& rOutputDir,
                                                                           const std::vector<std::string>& rOptions);

    /**
     * Get the pacing cycle length to use when simulating a cell.
     *
     * @param pCell  the cell to extract stimulus information from
     * @param defaultPeriod  default period if cell has no RegularStimulus
     */
    static double GetDefaultPeriod(boost::shared_ptr<AbstractCardiacCellInterface> pCell, double defaultPeriod);

    /**
     * Compare a simulation result against reference data, and compute a single error metric.
     * We use TODO!
     *
     * @param rSolution  the new simulation results
     * @param rModelName  the name of the model simulated, used to find the reference results in the project folder
     */
    static double GetError(const OdeSolution& rSolution, const std::string& rModelName);
};

#endif // CELLMODELUTILITIES_HPP_
