#pragma once

#include <blaze/Blaze.h>

class CartGrid
{
public:
	CartGrid(size_t nn);
	~CartGrid() {};

private:

	blaze::StaticVector<double, 3L> length_scales{ 1.0, 1.0, 1.0 };

	blaze::DynamicMatrix<int> grid_flags;
	blaze::DynamicMatrix<double> phi_matrix;
};