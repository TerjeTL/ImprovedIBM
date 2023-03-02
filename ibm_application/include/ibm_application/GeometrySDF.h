#pragma once

#include <Eigen/Dense>

//#include "ibm_application/SDLGraphics.h"
//#include "SDL.h"

class SDLGraphics;

enum class BoundaryCondition {
	Dirichlet,
	Neumann
};

class GeometrySDF
{
public:
	std::string name_id = "";

	GeometrySDF();
	GeometrySDF(double pos_x, double pos_y, double boundary_phi, bool inverse_sign = false);
	~GeometrySDF() = default;

	void SetPosition(double pos_x, double pos_y) { m_pos_x = pos_x; m_pos_y = pos_y; };
	void SetBoundaryCondition(BoundaryCondition boundary_condition)
	{
		m_boundary_condition = boundary_condition;
	}
	virtual double SignedDistanceFunction(double sample_x, double sample_y) const { return 0.0; };
	//virtual void RenderSDF(SDL_Renderer* renderer, SDLGraphics& graphics) {};

	Eigen::Vector2d GetPosition()
	{
		return Eigen::Vector2d{ m_pos_x, m_pos_y };
	}
	BoundaryCondition GetBoundaryCondition() const
	{
		return m_boundary_condition;
	}

	virtual Eigen::Vector2d GetNormal(double sample_x, double sample_y) { return Eigen::Vector2d{}; };
	virtual double GetBoundaryPhi() { return m_boundary_phi; };
	virtual double GetSize() = 0;

protected:
	BoundaryCondition m_boundary_condition = BoundaryCondition::Dirichlet;
	bool m_inverse_sign;
	double m_boundary_phi = 0.0;

	double m_pos_x = 0.0;
	double m_pos_y = 0.0;
};

class Circle2D_SDF : public GeometrySDF
{
public:
	Circle2D_SDF(double pos_x, double pos_y, double radius, double boundary_phi, bool inverse_sign = false)
		: GeometrySDF(pos_x, pos_y, boundary_phi, inverse_sign), m_radius(radius) {};
	~Circle2D_SDF() = default;


	virtual double SignedDistanceFunction(double sample_x, double sample_y) const override;
	//virtual void RenderSDF(SDL_Renderer* renderer, SDLGraphics& graphics) override;

	Eigen::Vector2d GetNormal(double sample_x, double sample_y) override;
	double GetSize()
	{
		return m_radius;
	}

private:
	double m_radius;
};