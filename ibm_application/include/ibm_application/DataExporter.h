#pragma once

#include "Solver.h"
#include "CartGrid.h"
#include "GeometrySDF.h"

#include <filesystem>
#include <iostream>
#include <Eigen/Core>
#include <highfive/H5Easy.hpp>

class DataExporter
{
public:
	enum class LoggingConfig
	{
		Transient,
		Steady
	};

	DataExporter(std::filesystem::path file_path, LoggingConfig logging_config)
		:	m_output_path(file_path), m_file(H5Easy::File(m_output_path.string(), H5Easy::File::Overwrite)), m_logging_type(logging_config)
	{

	}
	~DataExporter() = default;

	LoggingConfig GetLoggingConfig() const { return m_logging_type; }

	void SetDataRef(std::shared_ptr<std::map<size_t, Solution>> solutions)
	{
		m_solutions = solutions;
	}

	void GenerateHeaderInfos()
	{
		/*for (auto const& [mesh_level, solution] : *m_solutions)
		{
			std::string curr_dir = root_dir + "/mesh_" + std::to_string(mesh_level) + "/time_dict";

			double dt = solution.m_dt;
			size_t size = solution.m_iteration;
			std::vector<double> t(size);
			std::generate(t.begin(), t.end(), [n = 1, &dt]() mutable { return dt * n++; });

			H5Easy::dump(m_file, curr_dir, t);
		}*/
	}

	void WriteGeometry(std::string name, Circle2D_SDF geometry, double radius)
	{
		std::string curr_dir = "/geometry/" + name + "/data";

		// radius, x, y, bc_phi
		auto pos = geometry.GetPosition();
		int neumann_bc = geometry.GetBoundaryCondition() == BoundaryCondition::Neumann;
		std::vector<double> data = { radius, pos.x(), pos.y(), geometry.GetBoundaryPhi() };
		H5Easy::dump(m_file, curr_dir, data);

		curr_dir = "/geometry/" + name + "/bc";
		H5Easy::dump(m_file, curr_dir, neumann_bc);
	}

	void AppendCurrentState()
	{
		for (auto const& [mesh_level, solution] : *m_solutions)
		{
			if (solution.converged)
			{
				continue;
			}

			std::string curr_dir = root_dir + "/mesh_" + std::to_string(mesh_level) + time_dir + "/" + std::to_string(solution.m_time);

			// make a copy for now and store only active nodes
			auto mat = solution.m_mesh_grid->GetPhiMatrix();
			for (size_t i = 0; i < mat.rows(); i++)
			{
				for (size_t j = 0; j < mat.cols(); j++)
				{
					if (solution.m_mesh_grid->GetCellFlag(i, j) != 0)
					{
						mat(i, j) = 0;
					}
				}
			}

			H5Easy::dump(m_file, curr_dir, mat);

			// 2-norm
			/*curr_dir = root_dir + "/mesh_" + std::to_string(mesh_level) + "/euclidian_norm";
			double& data_value = 0.0;
			if (solution.m_iteration <= solution.m_iteration_level)
			{
				data_value = solution.euclidian_norm_init;
			}
			else
			{
				data_value = solution.euclidian_norm;
			}
			H5Easy::dump(m_file, curr_dir, solution.euclidian_norm, );*/
		}
	}

	void AppendSolutionData(const Solution& solution, int mesh_level, uint64_t logger_it)
	{
		// Time dict	
		std::string curr_dir = root_dir + "/mesh_" + std::to_string(mesh_level) + "/time_dict";
		H5Easy::dump(m_file, curr_dir, solution.m_time, { logger_it });


		curr_dir = root_dir + "/mesh_" + std::to_string(mesh_level) + time_dir + "/" + std::to_string(solution.m_time);
		// make a copy for now and store only active nodes
		auto mat = solution.m_mesh_grid->GetPhiMatrix();
		for (size_t i = 0; i < mat.rows(); i++)
		{
			for (size_t j = 0; j < mat.cols(); j++)
			{
				if (solution.m_mesh_grid->GetCellFlag(i, j) != 0)
				{
					mat(i, j) = 0;
				}
			}
		}

		H5Easy::dump(m_file, curr_dir, mat);

		// 2-norm
		curr_dir = root_dir + "/mesh_" + std::to_string(mesh_level) + "/euclidian_norm";
		double data_value = 0.0;
		if (solution.m_iteration == 1)
		{
			data_value = solution.euclidian_norm_init;
		}
		else
		{
			data_value = solution.euclidian_norm;
		}
		H5Easy::dump(m_file, curr_dir, data_value, { logger_it });
	
	}

	void WriteSteadyState(const Solution& solution, int mesh_level)
	{
		// End time	
		std::string curr_dir = root_dir + "/mesh_" + std::to_string(mesh_level) + "/steady_state/time";
		H5Easy::dump(m_file, curr_dir, solution.m_time);

		// Solution
		curr_dir = root_dir + "/mesh_" + std::to_string(mesh_level) + "/steady_state/solution";
		
		auto mat = solution.m_mesh_grid->GetPhiMatrix();
		for (size_t i = 0; i < mat.rows(); i++)
		{
			for (size_t j = 0; j < mat.cols(); j++)
			{
				if (solution.m_mesh_grid->GetCellFlag(i, j) != 0)
				{
					mat(i, j) = 0;
				}
			}
		}

		H5Easy::dump(m_file, curr_dir, mat);

		// Cell Flags
		curr_dir = root_dir + "/mesh_" + std::to_string(mesh_level) + "/steady_state/cell_flags";
		auto flags = solution.m_mesh_grid->GetGridFlags();
		H5Easy::dump(m_file, curr_dir, flags);

		// Track the first layer
		curr_dir = root_dir + "/mesh_" + std::to_string(mesh_level) + "/steady_state/boundary_values";
		
		Eigen::MatrixXd boundary_phi = Eigen::MatrixXd::Zero(mat.rows(), mat.cols());
		for (size_t i = 0; i < mat.rows(); i++)
		{
			for (size_t j = 0; j < mat.cols(); j++)
			{
				
				if (solution.m_mesh_grid->GetCellFlag(i, j) == 2)
				{
					auto east = solution.m_mesh_grid->GetCellFlag(i + 1, j);
					auto west = solution.m_mesh_grid->GetCellFlag(i - 1, j);
					auto north = solution.m_mesh_grid->GetCellFlag(i, j + 1);
					auto south = solution.m_mesh_grid->GetCellFlag(i, j - 1);

					if (east == 0)
					{
						boundary_phi(i + 1, j) = mat(i + 1, j);
					}
					if (west == 0)
					{
						boundary_phi(i - 1, j) = mat(i - 1, j);
					}
					if (north == 0)
					{
						boundary_phi(i, j + 1) = mat(i, j + 1);
					}
					if (south == 0)
					{
						boundary_phi(i, j - 1) = mat(i, j - 1);
					}
				}
			}
		}

		H5Easy::dump(m_file, curr_dir, boundary_phi);
	}

	void AppendMatrixData(std::string dir, const Eigen::MatrixXd& mat)
	{
		H5Easy::dump(m_file, dir, mat);
	}
private:
	LoggingConfig m_logging_type;
	std::shared_ptr<std::map<size_t, Solution>> m_solutions = nullptr;

	std::string root_dir = "/solutions";
	std::string time_dir = "/time_data";

	std::filesystem::path m_output_path;
	H5Easy::File m_file;
};