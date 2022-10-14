#pragma once

#include "Schemes.h"
#include "CartGrid.h"

#include <memory>

// A solution holds a grid and it's selected scheme
struct Solution
{
	double m_dt = 0.0;
	double m_von_neumann_num = 0.0;
	int m_iteration_level = 1;
	double m_time = 0.0;
	int m_iteration = 0;
	std::unique_ptr<FTCS_Scheme> m_scheme;
	std::shared_ptr<CartGrid> m_mesh_grid;
};