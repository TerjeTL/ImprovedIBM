#pragma once

#include "DataStructs.h"
#include "Schemes.h"
#include "CartGrid.h"
#include "RichardsonMethod.h"

#include <memory>
#include <chrono>

class DataExporter;

class Solver
{
public:
	Solver(double dt)
		: m_dt(dt)
	{
		m_time = m_start_time;
	};
	~Solver() {};

	void PerformStep(int steps = 1);
	void CheckConvergence(Solution& solution);
	void TaskStartPrintout(int task_iterations);
	void TaskFinishedPrintout();
	void SetDataExporter(std::shared_ptr<DataExporter> data_exporter) { m_data_export = data_exporter; };
	void SetRichardsonMethod(std::shared_ptr<RichardsonMethod> richardson_extrapolator) { m_richardson_extrapolator = richardson_extrapolator; };

	void AddSolution(size_t refinement_level, std::unique_ptr<FTCS_Scheme> scheme, std::shared_ptr<CartGrid> grid_mesh)
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

		(*m_solutions)[refinement_level] = Solution{ dt, von_neumann_num, iteration_level, m_time, 0, 0.0, false, std::move(scheme), grid_mesh };
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
	std::shared_ptr<RichardsonMethod> m_richardson_extrapolator = nullptr;
	std::shared_ptr<DataExporter> m_data_export = nullptr;

	bool m_converged = false;

	double m_alpha = 1.0;

	double m_start_time = 0.0;
	double m_end_time = 100.0;
	double m_time;
	double m_dt;
	double m_von_neumann_num = 0.0;
	int m_iterations = 0;

	double m_tolerance = 1.0e-22;
	double m_euclidian_norm_first_it = 0.0;

	std::chrono::duration<double> m_execution_time{ 0.0 };
};