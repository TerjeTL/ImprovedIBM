//#define SDL_MAIN_HANDLED
#include <SDL.h>

#include "data_viewer/SolutionModel.h"

#include "ibm_application/CartGrid.h"
#include "ibm_application/DataExporter.h"
#include "imgui.h"
#include "data_viewer/Camera.h"
#include "ibm_application/Solver.h"

class RichardsonExtrpGroup
{
public:
	RichardsonExtrpGroup(std::shared_ptr<Solver> solver) : m_solver(solver) {};

	void Update()
	{
		// Synchronize all solutions
		/*for (auto it = richardson_extrp_group.rbegin(); it != richardson_extrp_group.rend(); ++it)
		{
			const size_t size = it->first;
			m_solver->GetSolution(size)->CopyStateFromRefined();
		}*/

		auto solution = m_solver->GetSolution(richardson_extrp_group.begin()->first);
		solution->CopyStateFromRefined();

		// Update recursively starting from coarsest grid
		m_solver->GetSolution(richardson_extrp_group.begin()->first)->RecursiveUpdateFromThis();

		//Perform RE
		for (auto& [size, richardson_grid] : richardson_extrp_group)
		{
			if (size != richardson_extrp_group.rbegin()->first)
			{
				const Eigen::MatrixXd& curr_phi = m_solver->GetSolution(size)->m_mesh_grid->GetPhiMatrix();
				const Eigen::MatrixXd& fine_phi = m_solver->GetSolution((size*2-1))->m_mesh_grid->GetPhiMatrix();

				for (size_t i = 0; i < curr_phi.cols(); i++)
				{
					for (size_t j = 0; j < curr_phi.rows(); j++)
					{
						richardson_grid(j, i) = fine_phi(j * 2, i * 2) + (fine_phi(j * 2, i * 2) - curr_phi(j, i))/3.0;
					}
				}
			}
		}
	}

	void AddSolution(size_t grid_size)
	{
		if (m_solver->GetSolutions().find(grid_size) != m_solver->GetSolutions().end())
		{
			richardson_extrp_group[grid_size] = Eigen::MatrixXd::Zero(grid_size, grid_size);
		}
	}

	std::shared_ptr<Solver> m_solver;
	std::map<size_t, Eigen::MatrixXd> richardson_extrp_group;
	ImColor color = ImColor{ 1, 0, 0 };
private:
};

class DataViewer
{
public:
	DataViewer() : window_width( 800 ), window_height( 600 ) {}
	void DataViewerInitialize();
	~DataViewer() = default;

	void RunDataViewer();

	std::vector<RichardsonExtrpGroup> m_re_group;
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
	ImVector<size_t> table_selection;

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
