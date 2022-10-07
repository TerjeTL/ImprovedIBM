#pragma once

#include "Schemes.h"
#include "CartGrid.h"
#include "RichardsonMethod.h"

#include <memory>
#include <chrono>

class DataExporter;

// A solution holds a grid and it's selected scheme
struct Solution
{
	double m_dt = 0.0;
	std::unique_ptr<FTCS_Scheme> m_scheme;
	std::shared_ptr<CartGrid> m_mesh_grid;
};

class Solver
{
public:
	Solver(double dt, std::unique_ptr<FTCS_Scheme> scheme, std::shared_ptr<CartGrid> mesh)
		: m_dt(dt), m_selected_scheme(std::move(scheme)), m_grid_mesh(mesh)
	{
		m_time = m_start_time;
		auto h = m_grid_mesh->GetGridCellSize()(0);
		m_von_neumann_num = 2.0 * m_alpha * dt / (h*h);
	};
	~Solver() {};

	void PerformStep(int steps = 1);
	void CheckConvergence();
	void TaskStartPrintout(int task_iterations);
	void TaskFinishedPrintout();
	void SetDataExporter(std::shared_ptr<DataExporter> data_exporter) { m_data_export = data_exporter; };
	void SetRichardsonMethod(std::shared_ptr<RichardsonMethod> richardson_extrapolator) { m_richardson_extrapolator = richardson_extrapolator; };

	void AddSolution(size_t refinement_level, std::unique_ptr<FTCS_Scheme>, std::shared_ptr<CartGrid> grid)
	{
		auto dt = m_von_neumann_num * m_grid_mesh->GetGridCellSize()(0) / m_alpha
		m_solutions[refinement_level] = Solution{}
	}

	void AddMeshGrid(size_t refinement_level, std::shared_ptr<CartGrid> mesh_grid)
	{
		m_mesh_grids[refinement_level] = mesh_grid;
	}

	double GetCurrentTime() const
	{
		return m_time;
	}

private:
	std::unique_ptr<FTCS_Scheme> m_selected_scheme;
	std::map<size_t, Solution> m_solutions;
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

	double m_tolerance = 1.0e-8;
	double m_euclidian_norm_first_it = 0.0;

	std::chrono::duration<double> m_execution_time{ 0.0 };
};