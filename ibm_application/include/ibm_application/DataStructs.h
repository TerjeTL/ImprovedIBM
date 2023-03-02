#pragma once

#include "Schemes.h"
#include "CartGrid.h"

#include <iostream>
#include <memory>

// A solution holds a grid and it's selected scheme
struct Solution
{
	Solution() = default;

	Solution(double time, double dt, double von_neumann_num, std::unique_ptr<FTCS_Scheme> scheme, std::shared_ptr<CartGrid> mesh_grid)
		: m_time(time), m_dt(dt), m_von_neumann_num(von_neumann_num),
		m_scheme(std::move(scheme)), m_mesh_grid(std::move(mesh_grid))
	{
	
	}

	void Update()
	{
		m_scheme->Update(m_dt, m_von_neumann_num);
		m_iteration++;
		m_time += m_dt;
	}

	void TaskFinishedPrintout()
	{
		std::ios oldState(nullptr);
		oldState.copyfmt(std::cout);
		//--------------------------

		std::cout << std::boolalpha;

		std::cout << "\n======TASK FINISHED=========\n"
			<< "mesh size: " << m_mesh_grid->GetMeshSize().first << "x" << m_mesh_grid->GetMeshSize().second << "\n"
			<< "simulation time: " << m_time << "\n"
			<< "iterations: " << m_iteration << "\n\n";

		//--------------------------
		std::cout.copyfmt(oldState);
	}

	int m_iteration = 0;
	double m_time = 0.0;
	double m_dt = 0.0;
	double m_von_neumann_num = 0.0;

	std::shared_ptr<Solution> coarse_grid = nullptr;
	std::shared_ptr<Solution> fine_grid = nullptr;

	//double euclidian_norm = 0.0;
	//double euclidian_norm_init = 0.0;

	//int m_stop_iteration = -1;
	//bool converged = false;

	std::unique_ptr<FTCS_Scheme> m_scheme = nullptr;
	std::shared_ptr<CartGrid> m_mesh_grid = nullptr;
};