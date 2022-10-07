#pragma once

#include "Solver.h"
#include "CartGrid.h"

#include <filesystem>
#include <iostream>
#include <Eigen/Core>
#include <highfive/H5Easy.hpp>

class DataExporter
{
public:
	DataExporter(std::shared_ptr<Solver> solver)
		: m_file(H5Easy::File("export_data.h5", H5Easy::File::Overwrite)), m_solver(solver)
	{
		m_output_path = std::filesystem::current_path().parent_path().parent_path() / "scripts\\export_data.h5";
		m_file = H5Easy::File(m_output_path.string(), H5Easy::File::Overwrite);
	}
	~DataExporter() = default;


	void AppendCurrentState()
	{
		for (auto const& [mesh_level, solution] : m_solver->GetSolutions())
		{
			std::string curr_dir = root_dir + "/mesh_" + std::to_string(mesh_level) + time_dir + "/" + std::to_string(solution.m_time);
			H5Easy::dump(m_file, curr_dir, solution.m_mesh_grid->GetPhiMatrix());
		}
	}

	void AppendMatrixData(std::string dir, const Eigen::MatrixXd& mat)
	{
		H5Easy::dump(m_file, dir, mat);
	}
private:
	std::shared_ptr<Solver> m_solver = nullptr;

	std::string root_dir = "/solutions";
	std::string time_dir = "/time_data";

	std::filesystem::path m_output_path;
	H5Easy::File m_file;
};