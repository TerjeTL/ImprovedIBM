
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
		m_mesh_grids[refinement_level] = mesh_grid;
	}

	void ApplyExtrapolation()
	{

	}

private:
	std::map<size_t, std::shared_ptr<CartGrid>> m_mesh_grids;
};