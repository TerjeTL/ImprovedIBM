//#define SDL_MAIN_HANDLED
#pragma once
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
	RichardsonExtrpGroup(std::shared_ptr<Solver> solver, std::shared_ptr<Solution> fine_grid_solution) : m_solver(solver), m_fine_solution(fine_grid_solution)
	{
		size_t fine_grid_size = m_fine_solution->m_mesh_grid->GetMeshSize().first;

		size_t grid_size_coarse = (fine_grid_size + 1) / 2;
		double dt = fine_grid_solution->m_dt * 4.0;

		std::shared_ptr<CartGrid> coarse_grid_mesh{ std::make_shared<CartGrid>(grid_size_coarse) };
		auto h = coarse_grid_mesh->GetGridCellSize()(0);
		auto von_neumann_num = 2.0 * m_solver->GetThermalConductivity() * dt / (h * h);

		m_coarse_solution = std::make_shared<Solution>(m_fine_solution->m_time, dt, von_neumann_num, std::make_unique<FTCS_Scheme>(coarse_grid_mesh), std::move(coarse_grid_mesh));
		m_fine_solution->richardson_grid = m_coarse_solution;

		richardson_extrp = Eigen::MatrixXd::Zero(grid_size_coarse, grid_size_coarse);
	}

	void Update()
	{
		// Update recursively from fine grid
		m_fine_solution->RecursiveUpdateFromThis();

		//Perform RE
		const Eigen::MatrixXd& coarse_phi = m_coarse_solution->m_mesh_grid->GetPhiMatrix();
		const Eigen::MatrixXd& fine_phi = m_fine_solution->m_mesh_grid->GetPhiMatrix();

		for (size_t i = 0; i < coarse_phi.cols(); i++)
		{
			for (size_t j = 0; j < coarse_phi.rows(); j++)
			{
				richardson_extrp(j, i) = fine_phi(j * 2, i * 2) + (fine_phi(j * 2, i * 2) - coarse_phi(j, i))/3.0;
			}
		}
	}

	Eigen::MatrixXd GetRichardsonExtrapolation() const { return richardson_extrp; }

	std::shared_ptr<Solver> m_solver;

	std::shared_ptr<Solution> m_fine_solution;
	std::shared_ptr<Solution> m_coarse_solution;

	Eigen::MatrixXd richardson_extrp;
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

	std::map<size_t, RichardsonExtrpGroup> m_re_group;
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
