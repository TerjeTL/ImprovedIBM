#pragma once

#include <Eigen/Dense>

//#include "ibm_application/SDLGraphics.h"
#include "SDL.h"

class SDLGraphics;

class GeometrySDF
{
public:
	GeometrySDF();
	GeometrySDF(double pos_x, double pos_y, double boundary_phi);
	~GeometrySDF() = default;

	void SetPosition(double pos_x, double pos_y) { m_pos_x = pos_x; m_pos_y = pos_y; };
	virtual double SignedDistanceFunction(double sample_x, double sample_y) const { return 0.0; };
	virtual void RenderSDF(SDL_Renderer* renderer, SDLGraphics& graphics) {};
	virtual Eigen::Vector2d GetNormal(double sample_x, double sample_y) { return Eigen::Vector2d{}; };
	virtual double GetBoundaryPhi() { return m_boundary_phi; };

protected:
	double m_boundary_phi = 0.0;

	double m_pos_x = 0.0;
	double m_pos_y = 0.0;
};

class Circle2D_SDF : public GeometrySDF
{
public:
	Circle2D_SDF(double pos_x, double pos_y, double boundary_phi, double radius) : GeometrySDF(pos_x, pos_y, boundary_phi), m_radius(radius) {};
	~Circle2D_SDF() = default;


	virtual double SignedDistanceFunction(double sample_x, double sample_y) const override;
	virtual void RenderSDF(SDL_Renderer* renderer, SDLGraphics& graphics) override;

	Eigen::Vector2d GetNormal(double sample_x, double sample_y) override;

private:
	double m_radius;
};