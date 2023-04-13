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
		if (m_iteration == 0)
		{
			m_l2_norm_init = m_mesh_grid->GetPhiMatrix().squaredNorm() * m_mesh_grid->GetGridCellSize().x();
		}

		m_scheme->Update(m_dt, m_von_neumann_num);

		m_l2_norm_prev = m_l2_norm;
		m_l2_norm = m_mesh_grid->GetPhiMatrix().squaredNorm() * m_mesh_grid->GetGridCellSize().x();

		if (m_iteration == 0)
		{
			m_l2_norm_init = m_l2_norm / m_l2_norm_init;
			m_l2_norms.push_back(m_l2_norm_init);
		}
		else // need to skip calculation of residual when m_l2_norm_init = 0
		{
			m_l2_norms.push_back((m_l2_norm / m_l2_norm_prev) / m_l2_norm_init);
		}

		m_iteration++;
		m_time += m_dt;

		if (richardson_grid)
		{
			if (m_iteration % 4 == 0)
			{
				richardson_grid->CopyStateFromRefined();
				richardson_grid->Update();
			}
		}
	}

	void RecursiveUpdateFromThis()
	{
		Update();

		if (fine_grid)
		{
			fine_grid->RecursiveUpdate(this);
		}
		if (coarse_grid)
		{
			coarse_grid->RecursiveUpdate(this);
		}
		if (richardson_grid)
		{
			richardson_grid->RecursiveUpdate(this);
		}
	}

	void RecursiveUpdate(Solution* from_solution, int recursive_level = 1)
	{
		auto num_it = std::round(std::pow(4, recursive_level));

		if (coarse_grid.get() == from_solution) // func called from coarse grid
		{
			for (size_t i = 0; i < num_it; i++)
			{
				Update();
			}

			if (fine_grid)
			{
				fine_grid->RecursiveUpdate(this, ++recursive_level);
			}
		}
		else if (fine_grid.get() == from_solution) // func called from fine grid
		{
			fine_grid_it_delta++;
			if (fine_grid_it_delta >= num_it)
			{
				Update();
				fine_grid_it_delta = 0;
			}

			if (coarse_grid)
			{
				coarse_grid->RecursiveUpdate(this, ++recursive_level);
			}
		}
	}

	void CopyStateFromRefined()
	{
		if (fine_grid)
		{
			Eigen::MatrixXd& curr_grid_phi = m_mesh_grid->GetPhiMatrixRef();
			const Eigen::MatrixXd& fine_grid_phi = fine_grid->m_mesh_grid->GetPhiMatrixRef();

			for (size_t i = 0; i < curr_grid_phi.cols(); i++)
			{
				for (size_t j = 0; j < curr_grid_phi.rows(); j++)
				{
					curr_grid_phi(j, i) = fine_grid_phi(j * 2, i * 2);
				}
			}
		}
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
	std::shared_ptr<Solution> richardson_grid = nullptr;
	unsigned int fine_grid_it_delta = 0; // need to update this for every n number of fine grid iterations 

	double m_l2_norm_prev = 0.0;
	double m_l2_norm = 0.0;
	double m_l2_norm_init = 0.0;
	std::vector<double> m_l2_norms;

	//int m_stop_iteration = -1;
	//bool converged = false;

	std::unique_ptr<FTCS_Scheme> m_scheme = nullptr;
	std::shared_ptr<CartGrid> m_mesh_grid = nullptr;
};