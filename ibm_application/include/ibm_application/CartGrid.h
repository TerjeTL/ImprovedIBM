#pragma once

#include <memory>
#include <functional>
#include <blaze/Blaze.h>

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
	blaze::StaticVector<double, 2L> GetImagePoint(size_t i, size_t j) const
	{
		return image_points(i, j);
	}


	std::pair<int, int> GetMeshSize() const
	{
		return std::pair<int, int>{grid_flags.columns(), grid_flags.rows()};
	}

	std::unordered_map<std::size_t, std::shared_ptr<GeometrySDF>> GetImmersedBoundaries() const
	{
		return immersed_boundaries;
	}
	

private:
	std::unordered_map<std::size_t, std::shared_ptr<GeometrySDF>> immersed_boundaries;

	double h = 0.0;

	blaze::StaticVector<double, 3L> length_scales{ 1.0, 1.0, 1.0 };

	// Grid flags
	//##############################################################
	// Describe whether a cell is active, inactive or a ghost point
	// inactive:	0
	// active:		1
	// ghost point:	2
	blaze::DynamicMatrix<int> grid_flags;
	blaze::CompressedMatrix<blaze::StaticVector<double, 2L>> image_points;

	
	blaze::DynamicMatrix<double> phi_matrix;
};