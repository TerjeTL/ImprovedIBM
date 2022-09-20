
#pragma once

#include "ibm_application/GeometrySDF.h"
#include "ibm_application/SDLGraphics.h"
#include <cmath>

GeometrySDF::GeometrySDF(double pos_x, double pos_y, double boundary_phi) : m_pos_x(pos_x), m_pos_y(pos_y), m_boundary_phi(boundary_phi)
{

}

double Circle2D_SDF::SignedDistanceFunction(double sample_x, double sample_y) const
{
	double p_x = sample_x - m_pos_x;
	double p_y = sample_y - m_pos_y;

	return std::sqrt(p_x*p_x + p_y*p_y) - m_radius;
}

Eigen::Vector2d Circle2D_SDF::GetNormal(double sample_x, double sample_y)
{
    double dx = sample_x - m_pos_x;
    double dy = sample_y - m_pos_y;

    Eigen::Vector2d normal{dx, dy};
    normal.normalize();

    return normal;
}

void Circle2D_SDF::RenderSDF(SDL_Renderer* renderer, SDLGraphics& graphics)
{
    SDL_Color color = {255,255,255,255};
    GridPosF test = graphics.GetGridPosF(m_pos_x, m_pos_y);
    ScreenSpacePosF position = graphics.GetScreenSpacePosF(test);

    SDL_Color original;
    SDL_GetRenderDrawColor(renderer, &original.r, &original.g, &original.b, &original.a);
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);

    int render_radius = (int)(m_radius * (float)(graphics.GetGridExtent().first-1) * (float)graphics.GetGridCellSize());

    const int32_t diameter = (render_radius * 2);

    int32_t x = (render_radius - 1);
    int32_t y = 0;
    int32_t tx = 1;
    int32_t ty = 1;
    int32_t error = (tx - diameter);

    while (x >= y)
    {
        //  Each of the following renders an octant of the circle
        SDL_RenderDrawPoint(renderer, position.x + x, position.y - y);
        SDL_RenderDrawPoint(renderer, position.x + x, position.y + y);
        SDL_RenderDrawPoint(renderer, position.x - x, position.y - y);
        SDL_RenderDrawPoint(renderer, position.x - x, position.y + y);
        SDL_RenderDrawPoint(renderer, position.x + y, position.y - x);
        SDL_RenderDrawPoint(renderer, position.x + y, position.y + x);
        SDL_RenderDrawPoint(renderer, position.x - y, position.y - x);
        SDL_RenderDrawPoint(renderer, position.x - y, position.y + x);

        if (error <= 0)
        {
            ++y;
            error += ty;
            ty += 2;
        }

        if (error > 0)
        {
            --x;
            tx += 2;
            error += (tx - diameter);
        }
    }

    SDL_SetRenderDrawColor(renderer, original.r, original.g, original.b, original.a);
}