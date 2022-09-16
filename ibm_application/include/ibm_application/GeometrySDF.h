#pragma once

#include <blaze/Math.h>

//#include "ibm_application/SDLGraphics.h"
#include "SDL.h"

class SDLGraphics;

class GeometrySDF
{
public:
	GeometrySDF();
	GeometrySDF(double pos_x, double pos_y);
	~GeometrySDF() = default;

	void SetPosition(double pos_x, double pos_y) { m_pos_x = pos_x; m_pos_y = pos_y; };
	virtual double SignedDistanceFunction(double sample_x, double sample_y) const { return 0.0; };
	virtual void RenderSDF(SDL_Renderer* renderer, SDLGraphics& graphics) {};
	virtual blaze::StaticVector<double, 2L> GetNormal(double sample_x, double sample_y) { return blaze::StaticVector<double, 2L>{}; };

protected:

	double m_pos_x = 0.0;
	double m_pos_y = 0.0;
};

class Circle2D_SDF : public GeometrySDF
{
public:
	Circle2D_SDF(double pos_x, double pos_y, double radius) : GeometrySDF(pos_x, pos_y), m_radius(radius) {};
	~Circle2D_SDF() = default;


	virtual double SignedDistanceFunction(double sample_x, double sample_y) const override;
	virtual void RenderSDF(SDL_Renderer* renderer, SDLGraphics& graphics) override;

	blaze::StaticVector<double, 2L> GetNormal(double sample_x, double sample_y) override;

private:
	double m_radius;
};