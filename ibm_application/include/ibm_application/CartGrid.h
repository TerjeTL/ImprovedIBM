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

	void UpdateGrid();

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
	Eigen::Vector2d GetImagePoint(size_t i, size_t j) const
	{
		return image_points.at({ i, j });
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

	double BilinearInterpolation(Eigen::Vector2d world_coordinate);
	

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
};