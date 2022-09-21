#pragma once

#include "Schemes.h"
#include "CartGrid.h"

#include <memory>

class Solver
{
public:
	Solver(double dt, std::unique_ptr<FTCS_Scheme> scheme, std::shared_ptr<CartGrid> mesh)
		: m_dt(dt), m_selected_scheme(std::move(scheme)), m_grid_mesh(mesh)
	{
		m_time = m_start_time;
		m_cfl = m_alpha * dt / (m_grid_mesh->GetGridCellSize()(0));
	};
	~Solver() {};

	void PerformStep(int steps = 1);

private:
	std::unique_ptr<FTCS_Scheme> m_selected_scheme;
	std::shared_ptr<CartGrid> m_grid_mesh;

	double m_alpha = 1.0;

	double m_start_time = 0.0;
	double m_end_time = 10.0;
	double m_time;
	double m_dt;
	double m_cfl = 0.0;
};