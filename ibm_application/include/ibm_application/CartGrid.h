#pragma once

#include <memory>
#include <functional>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "ibm_application/GeometrySDF.h"

class CartGrid
{
public:
	CartGrid(size_t nn, double phi_0 = 0.0);
	~CartGrid() {};

	void InitializeField();
	void UpdateGrid();

	template <typename T> int sgn(T val) {
		return (T(0) < val) - (val < T(0));
	}

	void AddImmersedBoundary(std::string name, std::shared_ptr<GeometrySDF> geometry)
	{
		if (geometry) {

			immersed_boundaries[std::hash<std::string>{}(name)] = geometry;
		}
	}

	int GetCellFlag(size_t i, size_t j) const
	{
		return grid_flags(i, j);
	}
	
	double GetPhi(size_t i, size_t j) const
	{
		return phi_matrix(i, j);
	}

	Eigen::Vector2d GetImagePoint(size_t i, size_t j) const
	{
		return image_points.at({ i, j });
	}
	double GetImagePointPhi(size_t i, size_t j) const
	{
		return phi_image_point_matrix(i, j);
	}
	double GetBoundaryPhi(size_t i, size_t j) const
	{
		return boundary_phi(i, j);
	}

	std::pair<int, int> GetMeshSize() const
	{
		return std::pair<int, int>{grid_flags.cols(), grid_flags.rows()};
	}

	std::unordered_map<std::size_t, std::shared_ptr<GeometrySDF>> GetImmersedBoundaries() const
	{
		return immersed_boundaries;
	}

	Eigen::Vector2d GetGridCoordinate(Eigen::Vector2d world_coordinate)
	{
		return world_coordinate / h;
	}

	Eigen::Vector2d GetWorldCoordinate(Eigen::Vector2d grid_coordinate)
	{
		return grid_coordinate * h;
	}

	Eigen::Vector3d GetGridCellSize()
	{
		return length_scales * h;
	}

	double BilinearInterpolation(size_t i, size_t j);

	const Eigen::MatrixXd& GetPhiMatrix() const
	{
		return phi_matrix;
	}

	BoundaryCondition GetBoundaryCondition(size_t i, size_t j) const
	{
		return immersed_boundaries.at(ghost_point_parent_sdf(i, j))->GetBoundaryCondition();
	}

	Eigen::MatrixXd& GetPhiMatrixRef()
	{
		return phi_matrix;
	}

	Eigen::MatrixXi& GetGridFlags()
	{
		return grid_flags;
	}

	Eigen::MatrixXd& GetPhiImagePointMatrixRef()
	{
		return phi_image_point_matrix;
	}
	
	// Debugging-oriented functions
	std::array<Eigen::Vector2d, 4> GetBilinearInterpSelection(size_t i, size_t j)
	{
		return bilinear_interp_selection.at({ i, j });
	}
private:
	std::unordered_map<std::size_t, std::shared_ptr<GeometrySDF>> immersed_boundaries;

	double h = 0.0;

	Eigen::Vector3d length_scales{ 1.0, 1.0, 1.0 };

	// Grid flags
	//##############################################################
	// Describe whether a cell is active, inactive or a ghost point
	// inactive:	0
	// active:		1
	// ghost point:	2
	Eigen::MatrixXi grid_flags;

	std::map<std::pair<int, int>, Eigen::Vector2d> image_points;

	Eigen::MatrixXd phi_matrix;
	Eigen::MatrixXd phi_image_point_matrix;
	Eigen::MatrixXd boundary_phi;
	Eigen::MatrixX<size_t> ghost_point_parent_sdf;


	// Debugging-oriented variables
	std::map<std::pair<int, int>, std::array<Eigen::Vector2d, 4>> bilinear_interp_selection;
};