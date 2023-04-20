#include "ibm_application/DataExporter.h"
#include "data_viewer/DataViewer.h"

#include <cmath>

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

void DataExporter::WriteRichardsonExtrapolationData(const RichardsonExtrpGroup& richardson_extrp, size_t size)
{
	std::string base_dir = root_dir + "/richardson_extrp" + "/transient" + "/group_size_" + std::to_string(size) + "/" + std::to_string(richardson_extrp.m_coarse_solution->m_time);

	// Coarse solution
	std::string curr_dir = base_dir + "/coarse";

	const auto& coarse_solution = richardson_extrp.m_coarse_solution->m_mesh_grid;
	auto mat = coarse_solution->GetPhiMatrix();
	for (size_t j = 0; j < mat.rows(); j++)
	{
		for (size_t i = 0; i < mat.cols(); i++)
		{
			if (coarse_solution->GetCellFlag(i, j) != 0)
			{
				mat(j, i) = 0;
			}
		}
	}

	H5Easy::dump(m_file, curr_dir, mat);

	// Fine solution
	curr_dir = base_dir + "/fine";

	const auto& fine_solution = richardson_extrp.m_fine_solution->m_mesh_grid;
	mat = fine_solution->GetPhiMatrix();
	for (size_t j = 0; j < mat.rows(); j++)
	{
		for (size_t i = 0; i < mat.cols(); i++)
		{
			if (fine_solution->GetCellFlag(i, j) != 0)
			{
				mat(j, i) = 0;
			}
		}
	}

	H5Easy::dump(m_file, curr_dir, mat);

	// Richardson extrapolation
	curr_dir = base_dir + "/richardson_extrp";
	H5Easy::dump(m_file, curr_dir, richardson_extrp.richardson_extrp);

	transient_re_timelevels.push_back(std::to_string(richardson_extrp.m_coarse_solution->m_time));
}

template<class T>
std::vector<T>np_array_to_vec(py::array_t<T> py_array)
{
	return std::vector<T>(py_array.data(), py_array.data() + py_array.size());
}

void DataExporter::WriteAnalyticalTransientSolutions(size_t size, double r_outer)
{
	py::scoped_interpreter guard{};

	py::function analytical_solution =
		py::reinterpret_borrow<py::function>(   // cast from 'object' to 'function - use `borrow` (copy) or `steal` (move)
			py::module::import("bessel_stuff").attr("CalculateSolution")  // import method "min_rosen" from python "module"
			);

	for (const auto& time : transient_re_timelevels)
	{
		std::string dir = root_dir + "/richardson_extrp" + "/transient" + "/analytical_size_" + std::to_string(size) + "/" + time;

		//py::module_ analytical_solution = py::module_::import("bessel_stuff");
		py::array_t<double> result = analytical_solution(std::stod(time), m_solver->GetThermalConductivity(), r_outer);
		std::vector<double> phi = np_array_to_vec(result);

		auto solution_mesh = m_solver->GetSolution(size)->m_mesh_grid;
		Eigen::MatrixXi solution_mesh_grid = solution_mesh->GetGridFlags();

		Eigen::MatrixXd analytical = Eigen::MatrixXd::Zero(solution_mesh_grid.rows(), solution_mesh_grid.cols());
		for (int j = 0; j < solution_mesh_grid.rows(); ++j)
		{
			for (int i = 0; i < solution_mesh_grid.cols(); ++i)
			{
				auto coordinate = solution_mesh->GetWorldCoordinate({ i, j });

				if (solution_mesh_grid(j, i) == 0)
				{
					auto r = std::sqrt(std::pow(coordinate.x() - 0.5, 2) + std::pow(coordinate.y() - 0.5, 2));

					auto dr = r_outer / (phi.size() - 1);
					int idx_min = static_cast<int>(r / dr);

					double phi_interp = std::lerp(phi[idx_min], phi[idx_min + 1], r - idx_min * dr);

					analytical(j, i) = phi_interp * 2.0;
				}
			}
		}

		H5Easy::dump(m_file, dir, analytical);
	}
}