//#define SDL_MAIN_HANDLED
#pragma once
#include <SDL.h>

#include "data_viewer/SolutionModel.h"

#include "ibm_application/CartGrid.h"
#include "ibm_application/DataExporter.h"
#include "imgui.h"
#include "data_viewer/Camera.h"
#include "ibm_application/Solver.h"

enum class SimulationCase
{
	SteadyState,
	Unsteady
};

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
		m_coarse_solution->fine_grid = m_fine_solution;
		m_fine_solution->richardson_grid = m_coarse_solution;


		richardson_extrp = Eigen::MatrixXd::Zero(grid_size_coarse, grid_size_coarse);
	}

	void UpdateRichardsonExtrp()
	{
		// Update recursively from fine grid
		for (size_t i = 0; i < 4; i++)
		{
			m_fine_solution->RecursiveUpdateFromThis();
		}

		const Eigen::MatrixXd& fine_phi = m_fine_solution->m_mesh_grid->GetPhiMatrix();

		const Eigen::MatrixXd& coarse_phi = m_fine_solution->coarse_grid->m_mesh_grid->GetPhiMatrix();
		if (m_time_level_reset)
		{
			//coarse_phi = m_coarse_solution->m_mesh_grid->GetPhiMatrix();
		}

		for (size_t i = 0; i < coarse_phi.cols(); i++)
		{
			for (size_t j = 0; j < coarse_phi.rows(); j++)
			{
				if (m_fine_solution->m_mesh_grid->GetCellFlag(i*2, j*2) == 0)
				{
					richardson_extrp(j, i) = fine_phi(j * 2, i * 2) + (fine_phi(j * 2, i * 2) - coarse_phi(j, i)) / 3.0;
				}
			}
		}

		richardson_extrp_fine = Eigen::MatrixXd::Zero(fine_phi.rows(), fine_phi.cols());
		for (size_t i = 0; i < richardson_extrp_fine.cols(); i++)
		{
			for (size_t j = 0; j < richardson_extrp_fine.rows(); j++)
			{
				if (m_fine_solution->m_mesh_grid->GetCellFlag(i,j) == 0)
				{
					int i_coarse = i / 2;
					int j_coarse = j / 2;
					if ((i + 1) % 2 == 0 && j % 2 != 0 && m_fine_solution->m_mesh_grid->GetCellFlag(i + 1, j) == 0)
					{ // even i + 1, odd j
						double c1 = richardson_extrp(j_coarse, i_coarse) - fine_phi(j, i);
						double c2 = richardson_extrp(j_coarse, i_coarse + 1) - fine_phi(j, i + 2);
						double c = 1.0 / 2.0 * (c1 + c2);

						richardson_extrp_fine(j, i + 1) = fine_phi(j, i + 1) + c;
					}

					if (i % 2 != 0 && (j + 1) % 2 == 0 && m_fine_solution->m_mesh_grid->GetCellFlag(i, j+1) == 0)
					{ // even j + 1, odd i
						double c1 = richardson_extrp(j_coarse, i_coarse) - fine_phi(j, i);
						double c2 = richardson_extrp(j_coarse + 1, i_coarse) - fine_phi(j + 2, i);
						double c = 1.0 / 2.0 * (c1 + c2);

						richardson_extrp_fine(j + 1, i) = fine_phi(j + 1, i) + c;
					}

					if ((i + 1) % 2 == 0 && (j + 1) % 2 == 0 && m_fine_solution->m_mesh_grid->GetCellFlag(i+1, j+1) == 0)
					{ // even j + 1, odd i
						double c1 = richardson_extrp(j_coarse, i_coarse) - fine_phi(j, i);
						double c2 = richardson_extrp(j_coarse, i_coarse + 1) - fine_phi(j, i + 2);
						double c3 = richardson_extrp(j_coarse + 1, i_coarse) - fine_phi(j + 2, i);
						double c4 = richardson_extrp(j_coarse + 1, i_coarse + 1) - fine_phi(j + 2, i + 2);
						double c = 1.0 / 4.0 * (c1 + c2 + c3 + c4);

						richardson_extrp_fine(j + 1, i + 1) = fine_phi(j + 1, i + 1) + c;
					}

					if (i % 2 != 0 && j % 2 != 0)
					{
						richardson_extrp_fine(j, i) = richardson_extrp(j_coarse, i_coarse);
					}
				}
			}
		}
	}

	Eigen::MatrixXd GetRichardsonExtrapolation() const { return richardson_extrp; }

	std::shared_ptr<Solver> m_solver;

	std::shared_ptr<Solution> m_fine_solution;
	std::shared_ptr<Solution> m_coarse_solution;

	Eigen::MatrixXd richardson_extrp;
	Eigen::MatrixXd richardson_extrp_fine;
	ImColor color = ImColor{ 1, 0, 0 };

	bool m_time_level_reset = false;
private:

};

class DataViewer
{
public:
	DataViewer(SimulationCase case_setup) : selected_case(case_setup), window_width( 800 ), window_height( 600 ) {}
	void DataViewerInitialize();
	~DataViewer() = default;

	void RunDataViewer();

	std::map<size_t, RichardsonExtrpGroup> m_re_group;
	std::shared_ptr<Solver> m_solver{std::make_shared<Solver>()};
	DataExporter m_data_export{ std::filesystem::current_path().parent_path().parent_path() / "scripts/export_data.h5", DataExporter::LoggingConfig::Steady };
	std::vector<SolutionModel> models;
	std::vector<std::shared_ptr<GeometrySDF>> m_boundaries;
private:
	SimulationCase selected_case;

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
