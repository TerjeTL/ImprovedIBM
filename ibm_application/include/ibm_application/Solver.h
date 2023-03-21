#pragma once

#include "DataStructs.h"
#include "Schemes.h"
#include "CartGrid.h"
#include <memory>
#include <chrono>

class DataExporter;

class Solver
{
public:
	Solver()
	{
		m_time = m_start_time;
	};
	~Solver() {};

	void PerformStep(int steps = 1);
	void CheckConvergence(Solution& solution);
	//void TaskStartPrintout(int task_iterations);
	//void TaskFinishedPrintout();

	void AddSolution(double dt, size_t grid_size)
	{
		std::shared_ptr<CartGrid> grid_mesh{ std::make_shared<CartGrid>(grid_size) };

		auto h = grid_mesh->GetGridCellSize()(0);
		auto von_neumann_num = 2.0 * m_alpha * dt / (h * h);

		m_solutions[grid_size] = std::make_shared<Solution>(m_time, dt, von_neumann_num, std::make_unique<FTCS_Scheme>(grid_mesh), std::move(grid_mesh));
	}

	size_t AddGridDoubledSolution(const Solution& solution)
	{
		size_t grid_size = solution.m_mesh_grid->GetMeshSize().first; // using cols as size (assumption: square mesh)

		size_t grid_size_double = 2 * grid_size - 1;
		double dt = solution.m_dt / 4.0;

		if (m_solutions.find(grid_size_double) == m_solutions.end())
		{
			AddSolution(dt, grid_size_double);

			// im sure this wont end in tears
			m_solutions[grid_size]->fine_grid = m_solutions[grid_size_double];
			m_solutions[grid_size_double]->coarse_grid = m_solutions[grid_size];

			// return new size (usage optional)
			return grid_size_double;
		}
		else
		{
			return 0; // invalid operation
		}
	}

	double GetCurrentTime() const
	{
		return m_time;
	}

	const std::map<size_t, std::shared_ptr<Solution>>& GetSolutions() const
	{
		return m_solutions;
	}

	std::shared_ptr<Solution> GetSolution(size_t size_id)
	{
		return m_solutions.at(size_id);
	}

private:
	std::map<size_t, std::shared_ptr<Solution>> m_solutions = std::map<size_t, std::shared_ptr<Solution>>();

	double m_alpha = 1.0;

	double m_start_time = 0.0;
	double m_end_time = 100.0;
	double m_time = m_start_time;

	//std::chrono::duration<double> m_execution_time{ 0.0 };
	//std::chrono::duration<double> m_projected_time{ 0.0 };
};