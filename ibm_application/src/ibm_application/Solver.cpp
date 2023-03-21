#include "ibm_application/Solver.h"
#include "ibm_application/DataExporter.h"

#include <iostream>

void Solver::PerformStep(int steps)
{
	// Set up data to export
	//m_data_export->SetDataRef(m_solutions);
	//TaskStartPrintout(steps);

	//auto start_time = std::chrono::high_resolution_clock::now();

	//if (steps == -1)
	//{
	//	// treat as full simuation
	//	steps = static_cast<int>((m_end_time - m_start_time)/m_dt) + 1; // +1 to ensure we request to simulate the full length of the simulation
	//}

	for (int i = 0; i < steps; i++)
	{
		/*if (m_time+m_dt > m_end_time || m_converged)
		{
			auto end_time = std::chrono::high_resolution_clock::now();

			m_execution_time += std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);

			TaskFinishedPrintout();
			break;
		}*/

		//m_time += m_dt;
		//m_iterations += 1;
		
		for (auto& [mesh_size, solution] : m_solutions)
		{
			solution->Update();
		}
		
		/*if (m_iterations % 1 == 0)
		{
			m_current_progress = (float)m_iterations / (float)m_reference_iterations * 100.0f;

			auto current_time = std::chrono::high_resolution_clock::now();
			auto lapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(current_time - start_time);

			auto slope = lapsed_time.count() / m_current_progress;
			auto projection = (100.f - m_current_progress) * slope;

			std::cout << "\r" << "Running Simulation... " << (int)m_current_progress << "%   " << "ETA: " << (int)projection << "s              " << std::flush;
		}*/
	}

	//m_data_export->GenerateHeaderInfos();
}

void Solver::CheckConvergence(Solution& solution)
{
	/*solution.euclidian_norm = solution.m_scheme->GetEuclidianNorm() * static_cast<double>(solution.m_iteration_level * solution.m_iteration_level * solution.m_iteration_level);

	if (solution.euclidian_norm <= m_tolerance)
	{
		solution.converged = true;
	}*/
}

//void Solver::TaskStartPrintout(int task_iterations)
//{
//	std::ios oldState(nullptr);
//	oldState.copyfmt(std::cout);
//	//--------------------------
//
//	std::cout << std::boolalpha;
//
//	std::cout << "\n============================\n"
//		<< "         Task Start         \n"
//		<< "============================\n"
//		<< "dt: " << m_dt << "\n"
//		<< "r: " << m_von_neumann_num << "\n"
//		<< "Task Start: " << m_time << "\n"
//		<< "Task End: " << m_end_time << "\n"
//		<< "Task Iterations: " << task_iterations << "\n"
//		<< "Current Iteration: " << m_iterations << "\n"
//		<< "Solver Tolerance: " << m_tolerance << "\n\n";
//
//	//--------------------------
//	std::cout.copyfmt(oldState);
//}

//void Solver::TaskFinishedPrintout()
//{
//	std::ios oldState(nullptr);
//	oldState.copyfmt(std::cout);
//	//--------------------------
//	
//	std::cout << std::boolalpha;
//
//	std::cout << "\n\n============================\n"
//		<< "simulation time: " << m_time << "\n"
//		<< "converged: " << m_converged << "\n"
//		<< "iterations: " << m_iterations << "\n"
//		<< "execution time: " << m_execution_time.count() << "\n\n";
//
//	//--------------------------
//	std::cout.copyfmt(oldState);
//}