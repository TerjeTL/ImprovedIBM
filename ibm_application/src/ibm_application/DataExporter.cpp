#include "ibm_application/DataExporter.h"
#include "data_viewer/DataViewer.h"

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
}