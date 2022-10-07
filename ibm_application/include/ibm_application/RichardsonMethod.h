
#pragma once

#include "CartGrid.h"

#include <map>
#include <memory>

class RichardsonMethod
{
public:
	RichardsonMethod() {};
	~RichardsonMethod() {};

	void AddMeshGrid(size_t refinement_level, std::shared_ptr<CartGrid> mesh_grid)
	{
		if (refinement_level == 0)
		{
			phi_r = Eigen::MatrixXd::Zero(42,42);
		}

		m_mesh_grids[refinement_level] = mesh_grid;
	}

	void ApplyExtrapolation()
	{
		auto& coarse = m_mesh_grids.at(0)->GetPhiMatrixRef();
		auto& fine = m_mesh_grids.at(1)->GetPhiMatrixRef();

		for (size_t i = 0UL; i < coarse.rows(); ++i) {
			for (size_t j = 0UL; j < coarse.cols(); ++j) {
				if (m_mesh_grids.at(0)->GetCellFlag(i, j) == 0)
				{
					phi_r(i, j) = coarse(i, j) + (coarse(i, j) - fine(i * 2, j * 2)) / 3.0;
				}
			}
		}
	}

	const Eigen::MatrixXd& GetPhiMatrix() const
	{
		return phi_r;
	}

private:
	std::map<size_t, std::shared_ptr<CartGrid>> m_mesh_grids;
	Eigen::MatrixXd phi_r;
};