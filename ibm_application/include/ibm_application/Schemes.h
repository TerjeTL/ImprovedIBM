#pragma once
#include "CartGrid.h"

class FTCS_Scheme
{
public:
	FTCS_Scheme(std::shared_ptr<CartGrid> mesh_grid) : m_mesh_grid(mesh_grid)
	{
	
	};
	~FTCS_Scheme() {};

	void Update(double dt, double cfl);

	// hardcode boundary step
	void BoundaryCondition();

	double GetEuclidianNorm() const
	{
		return euclidian_norm;
	}

private:
	std::shared_ptr<CartGrid> m_mesh_grid;
	
	Eigen::MatrixXd phi_old;

	double euclidian_norm = 0.0;
};
