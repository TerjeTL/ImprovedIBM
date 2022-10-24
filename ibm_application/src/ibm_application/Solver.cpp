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

		int converged_solutions = 0;
		for (auto& [mesh_level, solution] : *m_solutions)
		{
			if (!solution.converged)
			{
				for (size_t i = 0; i < solution.m_iteration_level; i++)
				{
					solution.m_scheme->Update(solution.m_dt, solution.m_von_neumann_num);

					solution.m_time += solution.m_dt;
					solution.m_iteration++;

					if (solution.m_iteration > 1)
					{
						CheckConvergence(solution);
					}
					else if (solution.m_iteration == 1)
					{
						solution.euclidian_norm_init = solution.m_scheme->GetEuclidianNorm();
						m_data_export->AppendSolutionData(solution, mesh_level, 0);
					}

					if (m_data_export && solution.m_iteration % m_log_interval == 0)
					{
						m_data_export->AppendSolutionData(solution, mesh_level, solution.m_iteration/m_log_interval);
					}
				}
			}

			if (solution.converged)
			{
				converged_solutions++;
			}
		}

		if (converged_solutions == (*m_solutions).size())
		{
			m_converged = true;
		}

		if (m_richardson_extrapolator)
		{
			m_richardson_extrapolator->ApplyExtrapolation();
			std::string dir = "/solutions/mesh_r/time_data/" + std::to_string(m_time);
			m_data_export->AppendMatrixData(dir, m_richardson_extrapolator->GetPhiMatrix());
		}

		if (m_data_export)
		{
			//m_data_export->AppendCurrentState();
		}
	}

	m_data_export->GenerateHeaderInfos();
}

void Solver::CheckConvergence(Solution& solution)
{
	solution.euclidian_norm = solution.m_scheme->GetEuclidianNorm() * static_cast<double>(solution.m_iteration_level * solution.m_iteration_level * solution.m_iteration_level);

	if (solution.euclidian_norm <= m_tolerance)
	{
		solution.converged = true;
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