
#pragma once

#define SDL_MAIN_HANDLED
#include "SDL.h"

class SDLGraphics
{
public:
	SDLGraphics();
	void SDLGraphicsInitialize();
	~SDLGraphics() = default;

	void RunSDLGraphics();

private:
	SDL_Window* window;
	SDL_Renderer* renderer;

	int grid_cell_size = 25;
	int grid_width = 29;
	int grid_height = 23;

	// + 1 so that the last grid lines fit in the screen.
	int window_width = (grid_width * grid_cell_size) + 1;
	int window_height = (grid_height * grid_cell_size) + 1;

	// Place the grid cursor in the middle of the screen.
	SDL_Rect grid_cursor = {
		(grid_width - 1) / 2 * grid_cell_size,
		(grid_height - 1) / 2 * grid_cell_size,
		grid_cell_size,
		grid_cell_size
	};

	// The cursor ghost is a cursor that always shows in the cell below the
	// mouse cursor.
	SDL_Rect grid_cursor_ghost = { grid_cursor.x, grid_cursor.y, grid_cell_size,
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