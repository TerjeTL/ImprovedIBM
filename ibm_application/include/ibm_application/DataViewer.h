#pragma once

#include "data_viewer/MatrixObject.h"

#include <SDL.h>

#include "CartGrid.h"
#include "imgui.h"

class MatrixData
{
public:
	MatrixData(std::string name, std::pair<int, int> dim, MatrixObject matrix)
		: tag(name),
		size(dim),
		matrix_visualization(matrix)
	{};
	MatrixData(std::pair<int, int> dim, MatrixObject matrix)
		:
		size(dim),
		matrix_visualization(matrix)
	{
		tag = std::to_string(dim.first) + "x" + std::to_string(dim.second);
	};

	std::string tag;
	std::pair<int, int> size;

	// render object
	MatrixObject matrix_visualization;
};

class DataViewer
{
public:
	DataViewer() : window_width( 720 ), window_height( 480 ) {};
	void DataViewerInitialize();
	~DataViewer() = default;

	void RunDataViewer();


	std::vector<MatrixData> grids;
private:
	SDL_WindowFlags window_flags = (SDL_WindowFlags)(SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
	SDL_Window* window;
	SDL_GLContext gl_context;

	// Set window size based on grid extents
	int window_width;
	int window_height;

	SDL_bool quit = SDL_FALSE;
	SDL_bool mouse_active = SDL_FALSE;
	SDL_bool mouse_hover = SDL_FALSE;

	// Our state
	bool show_demo_window = true;
	bool show_another_window = false;
	ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

	unsigned int selected_mat = 0;
};