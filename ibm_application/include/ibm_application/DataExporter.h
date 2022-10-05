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
	DataExporter(std::shared_ptr<Solver> solver, std::shared_ptr<CartGrid> grid)
		: m_file(H5Easy::File("export_data.h5", H5Easy::File::Overwrite)),
		m_grid(grid), m_solver(solver)
	{}

	void AppendCurrentState()
	{
		std::string curr_dir = time_dir + "/" + std::to_string(m_solver->GetCurrentTime());

		H5Easy::dump(m_file, curr_dir, m_grid->GetPhiMatrix());
	}
private:
	std::shared_ptr<Solver> m_solver = nullptr;
	std::shared_ptr<CartGrid> m_grid = nullptr;

	std::string root_dir = "/solution";
	std::string time_dir = root_dir + "/time_data";

	std::filesystem::path output_path;
	H5Easy::File m_file;
};