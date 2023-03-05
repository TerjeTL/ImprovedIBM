//#define SDL_MAIN_HANDLED
#include <SDL.h>

#include "data_viewer/SolutionModel.h"

#include "ibm_application/CartGrid.h"
#include "ibm_application/DataExporter.h"
#include "imgui.h"
#include "data_viewer/Camera.h"
#include "ibm_application/Solver.h"

class DataViewer
{
public:
	DataViewer() : window_width( 720 ), window_height( 480 ) {}
	void DataViewerInitialize();
	~DataViewer() = default;

	void RunDataViewer();

	std::shared_ptr<Solver> m_solver{std::make_shared<Solver>()};
	DataExporter m_data_export{ std::filesystem::current_path().parent_path().parent_path() / "scripts/export_data.h5", DataExporter::LoggingConfig::Steady };
	std::vector<SolutionModel> models;
	std::vector<std::shared_ptr<GeometrySDF>> m_boundaries;
private:
	SDL_WindowFlags window_flags = (SDL_WindowFlags)(SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
	SDL_Window* window;
	SDL_GLContext gl_context;

	// Set window size based on grid extents
	int window_width;
	int window_height;

	uint16_t selected_mat = 0;

	bool m_run_simulation = false;
	int m_interations_remaining = 0;

	SDL_bool quit = SDL_FALSE;
	SDL_bool mouse_active = SDL_FALSE;
	SDL_bool mouse_hover = SDL_FALSE;

	// Our state
	Camera camera;

	bool show_demo_window = true;
	bool show_another_window = false;
	ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
};
