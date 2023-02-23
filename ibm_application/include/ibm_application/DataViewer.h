//#define SDL_MAIN_HANDLED
#include <SDL.h>

#include "CartGrid.h"
#include "imgui.h"

class DataViewer
{
public:
	DataViewer() : window_width( 720 ), window_height( 480 ) {};
	void DataViewerInitialize();
	~DataViewer() = default;

	void RunDataViewer();

	std::shared_ptr<CartGrid> grid = nullptr;
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
};