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
	Solver(double dt, int reference_iterations)
		: m_dt(dt), m_reference_iterations(reference_iterations)
	{
		m_time = m_start_time;
	};
	~Solver() {};

	void PerformStep(int steps = 1);
	void CheckConvergence(Solution& solution);
	void TaskStartPrintout(int task_iterations);
	void TaskFinishedPrintout();
	void SetDataExporter(std::shared_ptr<DataExporter> data_exporter) { m_data_export = data_exporter; };

	void AddSolution(size_t refinement_level, std::unique_ptr<FTCS_Scheme> scheme, std::shared_ptr<CartGrid> grid_mesh, int iterations = -1)
	{
		auto h = grid_mesh->GetGridCellSize()(0);
		int iteration_level = std::pow(2, refinement_level * 2);
		auto dt = m_dt / static_cast<double>(iteration_level);
		auto von_neumann_num = 2.0 * m_alpha * dt / (h * h);

		if (refinement_level == 0)
		{
			m_dt = dt;
			m_von_neumann_num = von_neumann_num;
		}

		(*m_solutions)[refinement_level] = Solution(m_time, dt, von_neumann_num, iteration_level, std::move(scheme), grid_mesh, iterations);
	}

	double GetCurrentTime() const
	{
		return m_time;
	}

	const std::map<size_t, Solution>& GetSolutions() const
	{
		return *m_solutions;
	}

private:
	std::shared_ptr <std::map<size_t, Solution>> m_solutions = std::make_shared<std::map<size_t, Solution>>();
	std::shared_ptr<DataExporter> m_data_export = nullptr;

	bool m_converged = false;

	double m_alpha = 1.0;

	double m_start_time = 0.0;
	double m_end_time = 100.0;
	double m_time;
	double m_dt;
	double m_von_neumann_num = 0.0;
	int m_iterations = 0;
	int m_reference_iterations = 0;

	int m_log_interval = 10;
	double m_tolerance = 1.0e-4;
	double m_euclidian_norm_first_it = 0.0;

	float m_current_progress = 0.0;

	std::chrono::duration<double> m_execution_time{ 0.0 };
	std::chrono::duration<double> m_projected_time{ 0.0 };
};