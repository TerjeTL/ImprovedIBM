#pragma once
#include "CartGrid.h"

class FTCS_Scheme
{
public:
	FTCS_Scheme() {};
	~FTCS_Scheme() {};

private:
	std::shared_ptr<CartGrid> m_mesh_grid;
};
