
#include "ibm_application/Schemes.h"

#include <omp.h>

void FTCS_Scheme::BoundaryCondition()
{
	auto grid_extents = m_mesh_grid->GetMeshSize();
	Eigen::MatrixXd& phi = m_mesh_grid->GetPhiMatrixRef();

	#pragma omp parallel for num_threads(8)
	for (int i = 0; i < grid_extents.first; i++)
	{
		for (int j = 0; j < grid_extents.second; j++)
		{
			// skip if node is not a ghost point
			if (m_mesh_grid->GetCellFlag(i, j) != 2)
			{
				continue;
			}

			auto& ip_ref = m_mesh_grid->GetPhiImagePointMatrixRef();
			ip_ref(i, j) = m_mesh_grid->BilinearInterpolation(i, j);

			auto image_pt_loc = m_mesh_grid->GetGridCoordinate(m_mesh_grid->GetImagePoint(i, j));
			Eigen::Vector2d ghost_pt_loc{ i, j };

			auto dl = std::abs((image_pt_loc - ghost_pt_loc).norm());
			
			switch (m_mesh_grid->GetBoundaryCondition(i, j))
			{
			case BoundaryCondition::Dirichlet:
			{
				// GP = 2*BI - IP
				phi(i, j) = 2 * m_mesh_grid->GetBoundaryPhi(i, j) - ip_ref(i, j);
				break;
			}
			case BoundaryCondition::Neumann:
			{
				// GP = IP - dl * d/dn(phi)|BI
				phi(i, j) = ip_ref(i, j) - dl * m_mesh_grid->GetBoundaryPhi(i, j);
				break;
			}
			default:
				break;
			}

			
		}
	}
}

void FTCS_Scheme::Update(double dt, double cfl)
{
	BoundaryCondition();

	auto grid_extents = m_mesh_grid->GetMeshSize();

	phi_old = m_mesh_grid->GetPhiMatrix();
	
	Eigen::MatrixXd& phi = m_mesh_grid->GetPhiMatrixRef();

	#pragma omp parallel for num_threads(8)
	for (int i = 1; i < grid_extents.first-1; i++)
	{
		for (int j = 1; j < grid_extents.second-1; j++)
		{
			// skip if node is a ghost point/inactive
			if (m_mesh_grid->GetCellFlag(i,j) != 0)
			{
				continue;
			}

			phi(i, j) = phi_old(i, j) + 0.5 * cfl * (phi_old(i + 1, j) - 2 * phi_old(i, j) + phi_old(i - 1, j) 
					+ phi_old(i, j + 1) - 2 * phi_old(i, j) + phi_old(i, j - 1));
		}
	}

	euclidian_norm = (phi - phi_old).squaredNorm();
}