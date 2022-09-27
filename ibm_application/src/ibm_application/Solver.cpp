#include "ibm_application/Solver.h"

#include <iostream>

void Solver::PerformStep(int steps)
{
	TaskStartPrintout(steps);

	auto start_time = std::chrono::high_resolution_clock::now();

	if (steps == -1)
	{
		// treat as full simuation
		steps = static_cast<int>((m_end_time - m_start_time)/m_dt) + 1; // +1 to ensure we request to simulate the full length of the simulation
	}

	for (size_t i = 0; i < steps; i++)
	{
		if (m_time+m_dt > m_end_time || m_converged)
		{
			auto end_time = std::chrono::high_resolution_clock::now();

			m_execution_time += std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);

			TaskFinishedPrintout();
			break;
		}

		m_selected_scheme->Update(m_dt, m_von_neumann_num);
		m_time += m_dt;
		m_iterations += 1;

		if (m_iterations > 1)
		{
			CheckConvergence();
		}
		else if (m_iterations == 1)
		{
			m_euclidian_norm_first_it = m_selected_scheme->GetEuclidianNorm();
		}
	}
}

void Solver::CheckConvergence()
{
	double euclidian_norm = m_selected_scheme->GetEuclidianNorm();

	if (euclidian_norm <= m_tolerance * m_euclidian_norm_first_it)
	{
		m_converged = true;
	}
}

void Solver::TaskStartPrintout(int task_iterations)
{
	std::ios oldState(nullptr);
	oldState.copyfmt(std::cout);
	//--------------------------

	std::cout << std::boolalpha;

	std::cout << "\n============================\n"
		<< "         Task Start         \n"
		<< "============================\n"
		<< "dt: " << m_dt << "\n"
		<< "r: " << m_von_neumann_num << "\n"
		<< "Task Start: " << m_time << "\n"
		<< "Task End: " << m_end_time << "\n"
		<< "Task Iterations: " << task_iterations << "\n"
		<< "Current Iteration: " << m_iterations << "\n"
		<< "Solver Tolerance: " << m_tolerance << "\n";

	//--------------------------
	std::cout.copyfmt(oldState);
}

void Solver::TaskFinishedPrintout()
{
	std::ios oldState(nullptr);
	oldState.copyfmt(std::cout);
	//--------------------------
	
	std::cout << std::boolalpha;

	std::cout << "============================\n"
		<< "simulation time: " << m_time << "\n"
		<< "converged: " << m_converged << "\n"
		<< "iterations: " << m_iterations << "\n"
		<< "execution time: " << m_execution_time.count() << "\n";

	//--------------------------
	std::cout.copyfmt(oldState);
}