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

#ifndef TESTODESOLVINGTIMES_HPP_
#define TESTODESOLVINGTIMES_HPP_

#include <boost/shared_ptr.hpp>
#include <cxxtest/TestSuite.h>

#include "Timer.hpp"

#include "CvodeAdaptor.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RegularStimulus.hpp"
#include "Shannon2004.hpp"
#include "Shannon2004Cvode.hpp"

class TestOdeSolvingTimes : public CxxTest::TestSuite
{
public:
    // I want to see how much slower CVODE would be in a tissue situation (PDE timestep 0.01ms).
    void TestShannonSolvingTimes() throw (Exception)
    {
        // Set up a default solver and a stimulus
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver());
        boost::shared_ptr<CvodeAdaptor> p_cvode_adaptor(new CvodeAdaptor());
        boost::shared_ptr<AbstractStimulusFunction> p_stimulus(new RegularStimulus(-25,5,1000,1));

        boost::shared_ptr<AbstractCardiacCell> shannon_euler(new CellShannon2004FromCellML(p_solver,p_stimulus));
        boost::shared_ptr<AbstractCardiacCell> shannon_cvode_adaptor(new CellShannon2004FromCellML(p_cvode_adaptor,p_stimulus));
        boost::shared_ptr<AbstractCvodeCell> shannon_cvode(new CellShannon2004FromCellMLCvode(p_solver,p_stimulus));

        double solution_time = 1000;
        double pde_time_step = 0.01;

        // A standard Forward Euler Solve.
        shannon_euler->SetTimestep(0.0025);

        std::cout << "Timings for 1000ms of solve:\n";
        Timer::Reset();
        for (double start_time = 0; start_time < solution_time; start_time+=pde_time_step)
        {
            shannon_euler->SolveAndUpdateState(start_time, start_time+pde_time_step);
        }
        Timer::Print("1. Forward-Euler");

        // A standard CVODE adaptor solve
        p_cvode_adaptor->SetForceReset(true);
        Timer::Reset();
        for (double start_time = 0; start_time < solution_time; start_time+=pde_time_step)
        {
            shannon_cvode_adaptor->SolveAndUpdateState(start_time, start_time+pde_time_step);
        }
        Timer::Print("2. CVODE adaptor (with resetting)");

        // A standard CVODE adaptor solve
        p_cvode_adaptor->SetForceReset(false);
        Timer::Reset();
        for (double start_time = 0; start_time < solution_time; start_time+=pde_time_step)
        {
            shannon_cvode_adaptor->SolveAndUpdateState(start_time, start_time+pde_time_step);
        }
        Timer::Print("3. CVODE adaptor (no resetting)");

        // A standard native CVODE solve
        shannon_cvode->SetForceReset(true);
        Timer::Reset();
        for (double start_time = 0; start_time < solution_time; start_time+=pde_time_step)
        {
            shannon_cvode->SolveAndUpdateState(start_time, start_time+pde_time_step);
        }
        Timer::Print("4. CVODE native (with resetting)");

        // Switch off the 'reset' on native CVODE so it saves its internal variables.
        shannon_cvode->SetForceReset(false);
        Timer::Reset();
        for (double start_time = 0; start_time < solution_time; start_time+=pde_time_step)
        {
            shannon_cvode->SolveAndUpdateState(start_time, start_time+pde_time_step);
        }
        Timer::Print("5. CVODE native (no resetting)");
     }
};

#endif // TESTODESOLVINGTIMES_HPP_
