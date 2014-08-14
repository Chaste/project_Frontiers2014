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

#include <iostream>
#include <boost/foreach.hpp>

#include "CellModelUtilities.hpp"

#include "Timer.hpp"

#include "CvodeAdaptor.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "RegularStimulus.hpp"

#include "shannon_wang_puglisi_weber_bers_2004.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Opt.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004CvodeOpt.hpp"

class TestOdeSolvingTimes : public CxxTest::TestSuite
{
public:
	void TestListingModels() throw (Exception)
	{
		std::vector<FileFinder> models(CellModelUtilities::GetListOfModels());
		std::cout << "Available models:" << std::endl;
		BOOST_FOREACH(FileFinder& r_model, models)
		{
		    std::cout << "  " << r_model.GetLeafNameNoExtension() << std::endl;
		}
		TS_ASSERT_LESS_THAN(20u, models.size());
	}

	/**
	 * For 1000 APs.
	 */
	void TestShannonSolvingTimes() throw (Exception)
	{
	    FileFinder model("projects/Frontiers2014/cellml/shannon_wang_puglisi_weber_bers_2004.cellml", RelativeTo::ChasteSourceRoot);
	    std::vector<std::string> options;
	    OutputFileHandler handler("TestOdeSolvingTimes_TestShannon_Euler");

		// Set up a default solver and a stimulus
		boost::shared_ptr<AbstractIvpOdeSolver> p_euler_solver(new EulerIvpOdeSolver());
		boost::shared_ptr<AbstractIvpOdeSolver> p_rk2_solver(new RungeKutta2IvpOdeSolver());
		boost::shared_ptr<AbstractIvpOdeSolver> p_rk4_solver(new RungeKutta4IvpOdeSolver());
		boost::shared_ptr<CvodeAdaptor> p_cvode_adaptor(new CvodeAdaptor());
        boost::shared_ptr<AbstractStimulusFunction> p_stimulus(new RegularStimulus(-25,5,1000,1));

        boost::shared_ptr<AbstractCardiacCellInterface> shannon_euler = CellModelUtilities::CreateCellModel(model, handler, options);
        shannon_euler->SetStimulusFunction(p_stimulus);

//		boost::shared_ptr<AbstractCardiacCell> shannon_euler(new Cellshannon_wang_puglisi_weber_bers_2004FromCellML(p_euler_solver,p_stimulus));
		boost::shared_ptr<AbstractCardiacCell> shannon_euler_opt(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLOpt(p_euler_solver,p_stimulus));

		boost::shared_ptr<AbstractCardiacCell> shannon_rk2(new Cellshannon_wang_puglisi_weber_bers_2004FromCellML(p_rk2_solver,p_stimulus));
		boost::shared_ptr<AbstractCardiacCell> shannon_rk2_opt(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLOpt(p_rk2_solver,p_stimulus));

		boost::shared_ptr<AbstractCardiacCell> shannon_rk4(new Cellshannon_wang_puglisi_weber_bers_2004FromCellML(p_rk4_solver,p_stimulus));
		boost::shared_ptr<AbstractCardiacCell> shannon_rk4_opt(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLOpt(p_rk4_solver,p_stimulus));

		boost::shared_ptr<AbstractCardiacCell> shannon_cvode_adaptor(new Cellshannon_wang_puglisi_weber_bers_2004FromCellML(p_cvode_adaptor,p_stimulus));
		boost::shared_ptr<AbstractCardiacCell> shannon_cvode_adaptor_opt(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLOpt(p_cvode_adaptor,p_stimulus));
        
		// Solver is ignored by native CVODE cells.
		boost::shared_ptr<AbstractCvodeCell> shannon_cvode(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_euler_solver,p_stimulus));
		boost::shared_ptr<AbstractCvodeCell> shannon_cvode_opt(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvodeOpt(p_euler_solver,p_stimulus));

		double start_time = 0;
		double end_time = 10*1000;
		std::cout << "Timings for " << end_time/1000.0 << "s of 1Hz pacing:\n";

		///\todo IMPORTANT
		// figure out how to set the timesteps for the fixed methods 'fairly'
		// at present they are set so that all of them are just about stable.
		// But we probably want to set them for equivalent levels of accuracy.

		// A standard Forward Euler Solve.
		shannon_euler->SetTimestep(0.0025);
		Timer::Reset();
		shannon_euler->SolveAndUpdateState(start_time, end_time);
		Timer::Print("1. Forward-Euler");

		// An Opt Forward Euler Solve.
		shannon_euler_opt->SetTimestep(0.0025);
		// Set up lookup tables
		shannon_euler_opt->SolveAndUpdateState(-0.0025, 0);

		Timer::Reset();
		shannon_euler_opt->SolveAndUpdateState(start_time, end_time);
		Timer::Print("2. Forward-Euler Opt");

		// A standard RK2 Solve.
		shannon_rk2->SetTimestep(0.01);
		Timer::Reset();
		shannon_rk2->SolveAndUpdateState(start_time, end_time);
		Timer::Print("3. Runge-Kutta 2nd order");

		// An Opt RK2 Solve.
		shannon_rk2_opt->SetTimestep(0.01);
		Timer::Reset();
		shannon_rk2_opt->SolveAndUpdateState(start_time, end_time);
		Timer::Print("4. Runge-Kutta 2nd order Opt");

		// A standard RK4 Solve.
		shannon_rk4->SetTimestep(0.025);
		Timer::Reset();
		shannon_rk4->SolveAndUpdateState(start_time, end_time);
		Timer::Print("5. Runge-Kutta 4th order ");

		// An Opt RK4 Solve.
		shannon_rk4_opt->SetTimestep(0.025);
		Timer::Reset();
		shannon_rk4_opt->SolveAndUpdateState(start_time, end_time);
		Timer::Print("6. Runge-Kutta 4th order Opt");


		p_cvode_adaptor->SetMaxSteps(1e6);
		p_cvode_adaptor->SetTolerances(1e-5, 1e-7); // Match defaults for native cell
		shannon_cvode_adaptor->SetTimestep(boost::static_pointer_cast<RegularStimulus>(p_stimulus)->GetDuration());
		Timer::Reset();
		shannon_cvode_adaptor->SolveAndUpdateState(start_time, end_time);
		Timer::Print("7. CVODE Adaptor (no resetting)");

		shannon_cvode_adaptor_opt->SetTimestep(boost::static_pointer_cast<RegularStimulus>(p_stimulus)->GetDuration());
		Timer::Reset();
		shannon_cvode_adaptor_opt->SolveAndUpdateState(start_time, end_time);
		Timer::Print("8. CVODE Adaptor Opt (no resetting)");

		// A standard native CVODE solve
		shannon_cvode->SetMaxSteps(1e6);
		shannon_cvode->SetMaxTimestep(boost::static_pointer_cast<RegularStimulus>(p_stimulus)->GetDuration());

		TS_ASSERT_EQUALS(shannon_cvode->GetUseAnalyticJacobian(), true);
		shannon_cvode->ForceUseOfNumericalJacobian(true);
		TS_ASSERT_EQUALS(shannon_cvode->GetUseAnalyticJacobian(), false);

		Timer::Reset();
		shannon_cvode->SolveAndUpdateState(start_time, end_time);
		Timer::Print("9. CVODE Numerical Jacobian (native, no resetting)");

		shannon_cvode->ForceUseOfNumericalJacobian(false);
		TS_ASSERT_EQUALS(shannon_cvode->GetUseAnalyticJacobian(), true);
		shannon_cvode->ResetToInitialConditions();

		Timer::Reset();
		shannon_cvode->SolveAndUpdateState(start_time, end_time);
		Timer::Print("10. CVODE Analytic Jacobian (native, no resetting)");

		// A standard native CVODE solve
		shannon_cvode_opt->SetMaxSteps(1e6);
		shannon_cvode_opt->SetMaxTimestep(boost::static_pointer_cast<RegularStimulus>(p_stimulus)->GetDuration());

		TS_ASSERT_EQUALS(shannon_cvode_opt->GetUseAnalyticJacobian(), true);
		shannon_cvode_opt->ForceUseOfNumericalJacobian(true);
		TS_ASSERT_EQUALS(shannon_cvode_opt->GetUseAnalyticJacobian(), false);

		// IMPORTANT
		// An initial call to solve on Opt cells with CVODE seems to take about 0.9 seconds.
		// It might be setting up the lookup tables for later use, so
		// not really fair to time this, but we should make sure users know that they shouldn't
		// be setting up a new model each time with Opt cells.
		//
		// (doesn't seem to be the case for other sorts of cells - are there lookup tables for
		// the analytic jacobian entries that get calculated too?)
		shannon_cvode_opt->SolveAndUpdateState(-0.01, 0);

		Timer::Reset();
		shannon_cvode_opt->SolveAndUpdateState(start_time, end_time);
		Timer::Print("11. CVODE Opt Numerical Jacobian (native, no resetting)");

		shannon_cvode_opt->ForceUseOfNumericalJacobian(false);
		TS_ASSERT_EQUALS(shannon_cvode_opt->GetUseAnalyticJacobian(), true);
		shannon_cvode_opt->ResetToInitialConditions();

		Timer::Reset();
		shannon_cvode_opt->SolveAndUpdateState(start_time, end_time);
		Timer::Print("12. CVODE Opt Analytic Jacobian (native, no resetting)");
	}

	// I wanted to see how much slower CVODE would be in a tissue situation (pretend there is a PDE timestep 0.01ms).
	// result - it would be faster than Forward Euler!
	void TestShannonSolvingTimesForTissue() throw (Exception)
    {
		// Set up a default solver and a stimulus
		boost::shared_ptr<AbstractIvpOdeSolver> p_euler_solver(new EulerIvpOdeSolver());
		boost::shared_ptr<AbstractStimulusFunction> p_stimulus(new RegularStimulus(-25,5,1000,1));

		boost::shared_ptr<AbstractCardiacCell> shannon_euler(new Cellshannon_wang_puglisi_weber_bers_2004FromCellML(p_euler_solver,p_stimulus));
		boost::shared_ptr<AbstractCvodeCell> shannon_cvode(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_euler_solver,p_stimulus));

		double solution_time = 1000;
		double pde_time_step = 0.01;

		// A standard Forward Euler Solve.
		shannon_euler->SetTimestep(0.0025);

		for (unsigned i=0; i<2; i++)
		{
			if (i==0)
			{
				pde_time_step = 0.01;
			}
			else
			{
				pde_time_step = 0.1;
			}

			std::cout << "\nTimings for 1000ms of 'tissue' solve with pde_time_step = " << pde_time_step << std::endl;
			Timer::Reset();
			for (double start_time = 0; start_time < solution_time; start_time+=pde_time_step)
			{
				shannon_euler->SolveAndUpdateState(start_time, start_time+pde_time_step);
			}
			Timer::Print("1. Forward-Euler");

			// A standard native CVODE solve
			shannon_cvode->SetForceReset(true);
			shannon_cvode->SetMaxTimestep(boost::static_pointer_cast<RegularStimulus>(p_stimulus)->GetDuration());
			Timer::Reset();
			for (double start_time = 0; start_time < solution_time; start_time+=pde_time_step)
			{
				shannon_cvode->SolveAndUpdateState(start_time, start_time+pde_time_step);
			}
			Timer::Print("2. CVODE native (with resetting)");

			// Switch off the 'reset' on native CVODE so it saves its internal variables.
			shannon_cvode->SetForceReset(false);
			shannon_cvode->ResetToInitialConditions();
			Timer::Reset();
			for (double start_time = 0; start_time < solution_time; start_time+=pde_time_step)
			{
				shannon_cvode->SolveAndUpdateState(start_time, start_time+pde_time_step);
			}
			Timer::Print("3. CVODE native (no resetting)");
		}
    }

private:
	/* Utility methods used by the tests above go here.
	 */

//	/**
//	 * Find out how long it takes to simulate the given model.
//	 * The cell will be reset to initial conditions prior to simulation.
//	 * We assume the cell has a regular square wave stimulus defined.
//	 *
//	 * Note that we don't check the results are sensible, or do a pre-simulation to avoid counting lookup tables setup.
//	 *
//	 * @param pCell  the cell model to simulate, with solver attached
//	 * @param numPaces  the number of simulated paces to time
//	 * @return  elapsed wall clock time, in seconds
//	 */
//	double TimeSimulation(boost::shared_ptr<AbstractCardiacCellInterface> pCell,
//	                      unsigned numPaces)
//	{
//	    pCell->ResetToInitialConditions();
//	    double period = CellModelUtilities::GetDefaultPeriod(pCell, 1000.0);
//	    Timer::Reset();
//	    pCell->SolveAndUpdateState(0.0, numPaces * period);
//	    return Timer::GetElapsedTime();
//	}
};

#endif // TESTODESOLVINGTIMES_HPP_
