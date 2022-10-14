#include "ibm_application/Solver.h"
#include "ibm_application/DataExporter.h"

#include <iostream>

void Solver::PerformStep(int steps)
{
	// Set up data to export
	m_data_export->SetDataRef(m_solutions);
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

		m_time += m_dt;
		m_iterations += 1;

		for (auto& [mesh_level, solution] : *m_solutions)
		{
			for (size_t i = 0; i < solution.m_iteration_level; i++)
			{
				solution.m_scheme->Update(solution.m_dt, solution.m_von_neumann_num);
				solution.m_time += solution.m_dt;
				solution.m_iteration++;
			}
		}

		if (m_richardson_extrapolator)
		{
			m_richardson_extrapolator->ApplyExtrapolation();
			std::string dir = "/solutions/mesh_r/time_data/" + std::to_string(m_time);
			m_data_export->AppendMatrixData(dir, m_richardson_extrapolator->GetPhiMatrix());
		}
		

		if (m_iterations > 1)
		{
			CheckConvergence();
		}
		else if (m_iterations == 1)
		{
			m_euclidian_norm_first_it = m_solutions->at(1).m_scheme->GetEuclidianNorm();
		}

		if (m_data_export)
		{
			m_data_export->AppendCurrentState();
		}
	}

	m_data_export->GenerateHeaderInfos();
}

void Solver::CheckConvergence()
{
	double euclidian_norm = m_solutions->at(1).m_scheme->GetEuclidianNorm();

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
		<< "execution time: " << m_execution_time.count() << "\n\n";

	//--------------------------
	std::cout.copyfmt(oldState);
}