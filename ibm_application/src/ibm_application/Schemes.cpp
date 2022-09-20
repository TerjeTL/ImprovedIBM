
#include "ibm_application/Schemes.h"

void FTCS_Scheme::MakeIteration(double dt, double cfl)
{
	auto grid_extents = m_mesh_grid->GetMeshSize();

	phi_old = m_mesh_grid->GetPhiMatrix();
	
	Eigen::MatrixXd& phi = m_mesh_grid->GetPhiMatrixRef();

	for (size_t i = 0; i < grid_extents.first; i++)
	{
		for (size_t j = 0; j < grid_extents.second; j++)
		{
			// skip is node is a ghost point/inactive
			if (m_mesh_grid->GetCellFlag(i,j) != 0)
			{
				continue;
			}

			phi(i, j) = phi_old(i, j) + cfl * (phi_old(i + 1, j) - 2 * phi_old(i, j) + phi_old(i - 1, j) 
					+ phi_old(i, j + 1) - 2 * phi_old(i, j) + phi_old(i, j - 1));
		}
	}
}