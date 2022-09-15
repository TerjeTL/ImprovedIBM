
#pragma once

#define SDL_MAIN_HANDLED
#include "SDL.h"
#include "ibm_application/GeometrySDF.h"

#include <iostream>

struct GridPos
{
	int x, y;
};

struct GridPosF
{
	double x, y;
};

struct ScreenSpacePos
{
	int x, y;
};

struct ScreenSpacePosF
{
	double x, y;
};

class SDLGraphics
{
public:
	SDLGraphics();
	void SDLGraphicsInitialize();
	~SDLGraphics() = default;

	void RunSDLGraphics();

	ScreenSpacePos GetScreenSpacePos(GridPos grid_location);
	ScreenSpacePosF GetScreenSpacePosF(GridPosF grid_location);
	GridPosF GetGridPosF(double x_pos_norm, double y_pos_norm);
	int GetGridCellSize() { return grid_cell_size; };

private:
	SDL_Window* window;
	SDL_Renderer* renderer;

	Circle2D_SDF immersed_boundary{ 15., 15., 6. };

	void RenderCircle(SDL_Renderer* renderer, int32_t centreX, int32_t centreY, int32_t radius, SDL_Color color = { 255, 255, 255, 255 });
	void RenderFillCircle(SDL_Renderer* renderer, int32_t centreX, int32_t centreY, int32_t radius, SDL_Color color = { 255, 255, 255, 255 });

	int grid_cell_size = 25;
	int zoom_gain = 0;
	int zoom_level = 0;
	int grid_width = 31;
	int grid_height = 31;

	// Set window size based on grid extents
	int window_width = (grid_width * grid_cell_size) + 2;
	int window_height = (grid_height * grid_cell_size) + 2;

	// Place grid cursor in the centre of the window
	SDL_FRect grid_cursor = {
		((float)grid_width - 1.f) / 2.f * (float)grid_cell_size,
		((float)grid_height - 1.f) / 2.f * (float)grid_cell_size,
		(float)grid_cell_size,
		(float)grid_cell_size
	};
	
	// Select starting grid position at/near the centre
	GridPos grid_position = {
		(grid_width - 1) / 2,
		(grid_height - 1) / 2
	};

	// The cursor ghost is a cursor that always shows in the cell below the
	// mouse cursor.
	SDL_Rect grid_cursor_ghost = { 0, 0, grid_cell_size,
								  grid_cell_size };

	// Dark theme.
	SDL_Color grid_background = { 22, 22, 22, 255 }; // Barely Black
	SDL_Color grid_line_color = { 44, 44, 44, 255 }; // Dark grey
	SDL_Color grid_cursor_ghost_color = { 44, 44, 44, 255 };
	SDL_Color grid_cursor_color = { 255, 255, 255, 255 }; // White

	SDL_bool quit = SDL_FALSE;
	SDL_bool mouse_active = SDL_FALSE;
	SDL_bool mouse_hover = SDL_FALSE;
};