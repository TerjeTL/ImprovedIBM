
#include "ibm_application/Schemes.h"
#include <iostream>

#include <omp.h>

#define MT_ON
//#define DEBUG_TERMINAL_WLSQ

void FTCS_Scheme::BoundaryCondition()
{
	// Run WLSQ
	m_mesh_grid->WeightedLeastSquaresMethod();

	auto grid_extents = m_mesh_grid->GetMeshSize();
	Eigen::MatrixXd& phi = m_mesh_grid->GetPhiMatrixRef();

	#ifdef MT_ON
		#pragma omp parallel for num_threads(4)
	#endif
	for (int j = 0; j < grid_extents.first; j++)
	{
		for (int i = 0; i < grid_extents.second; i++)
		{
			// skip if node is not a ghost point
			if (m_mesh_grid->GetCellFlag(i, j) != 2)
			{
				continue;
			}

			//auto& ip_ref = m_mesh_grid->GetPhiImagePointMatrixRef();
			//ip_ref(j, i) = m_mesh_grid->BilinearInterpolation(i, j);

			//auto image_pt_loc_wrld = m_mesh_grid->GetImagePoint(i, j);
			//auto image_pt_loc_grid = m_mesh_grid->GetGridCoordinate(image_pt_loc_wrld);
			//Eigen::Vector2d ghost_pt_loc_grid{ i, j };
			//auto ghost_pt_loc_wrld = m_mesh_grid->GetWorldCoordinate(ghost_pt_loc_grid);

			//auto dr = image_pt_loc_wrld - ghost_pt_loc_wrld;

			//double dl = std::abs(dr.norm());

			phi(j,i) = m_mesh_grid->GetPhiWLSQ(i, j);

			//switch (m_mesh_grid->GetBoundaryCondition(i, j))
			//{
			//case BoundaryCondition::Dirichlet:
			//{
			//	// WLSQ Method
			//	phi(j, i) = m_mesh_grid->GetPhiWLSQ(i, j);

			//	// Image Point Method
			//	// GP = IP + (BI - IP)*len_factor
			//	//phi(j, i) = ip_ref(j, i) + (m_mesh_grid->GetBoundaryPhi(i, j) - ip_ref(j, i)) * m_mesh_grid->m_ip_stencil_length_factor;
			//	break;
			//}
			//case BoundaryCondition::Neumann:
			//{
			//	// WLSQ Method
			//	phi(j, i) = m_mesh_grid->GetPhiWLSQ(i, j);

			//	// Image Point Method
			//	// GP = IP - dl * d/dn(phi)|BI
			//	//phi(j, i) = ip_ref(j, i) - dl * m_mesh_grid->GetBoundaryPhi(i, j);
			//	break;
			//}
			//default:
			//	break;
			//}

		}
	}

#ifdef DEBUG_TERMINAL_WLSQ
	std::cout << "PHI MATRIX\n" << m_mesh_grid->GetPhiMatrix() << "\n\n";
#endif
}

void FTCS_Scheme::Update(double dt, double r)
{
	BoundaryCondition();

	auto grid_extents = m_mesh_grid->GetMeshSize();

	phi_old = m_mesh_grid->GetPhiMatrix();
	
	Eigen::MatrixXd& phi = m_mesh_grid->GetPhiMatrixRef();

	//std::cout << phi << "\n\n";

	#ifdef MT_ON
		#pragma omp parallel for num_threads(4)
	#endif
	for (int j = 1; j < grid_extents.first-1; j++)
	{
		for (int i = 1; i < grid_extents.second-1; i++)
		{
			// skip if node is a ghost point/inactive
			if (m_mesh_grid->GetCellFlag(i,j) != 0)
			{
				continue;
			}

			phi(j, i) = phi_old(j, i) + 1.0 * r * (phi_old(j, i + 1) - 2 * phi_old(j, i) + phi_old(j, i - 1) 
					+ phi_old(j + 1, i) - 2 * phi_old(j, i) + phi_old(j - 1, i));
		}
	}

	euclidian_norm = (phi - phi_old).squaredNorm();
}