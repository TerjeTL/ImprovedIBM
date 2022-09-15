#pragma once

//#include "ibm_application/SDLGraphics.h"
#include "SDL.h"

class SDLGraphics;

class GeometrySDF
{
public:
	GeometrySDF(double pos_x, double pos_y);
	~GeometrySDF() = default;

	void SetPosition(double pos_x, double pos_y) { m_pos_x = pos_x; m_pos_y = pos_y; };
	virtual double SignedDistanceFunction(double sample_x, double sample_y) = 0;
	virtual void RenderSDF(SDL_Renderer* renderer, SDLGraphics& graphics) = 0;

protected:

	double m_pos_x;
	double m_pos_y;
};

class Circle2D_SDF : public GeometrySDF
{
public:
	Circle2D_SDF(double pos_x, double pos_y, double radius) : GeometrySDF(pos_x, pos_y), m_radius(radius) {};
	~Circle2D_SDF() = default;


	virtual double SignedDistanceFunction(double sample_x, double sample_y) override;
	virtual void RenderSDF(SDL_Renderer* renderer, SDLGraphics& graphics) override;
private:

	double m_radius;
};