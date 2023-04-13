#pragma once

#include "data_viewer/DataViewer.h"
#include "data_viewer/Shader.h"
#include "data_viewer/Camera.h"
#include "data_viewer/SolutionModel.h"

#include "glad/glad.h"
#include "imgui.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"
#include <stdio.h>
#include <SDL.h>
#include <SDL_opengl.h>
#include <vector>
#include "glm.hpp"
#include "imgui_internal.h"
#include "implot_internal.h"

#include "Eigen/Eigen"

#include "data_viewer/Mesh.h"

void GLAPIENTRY
MessageCallback(GLenum source,
    GLenum type,
    GLuint id,
    GLenum severity,
    GLsizei length,
    const GLchar* message,
    const void* userParam)
{
    if (severity == GL_DEBUG_SEVERITY_HIGH)
	{
        fprintf(stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
            (type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""),
            type, severity, message);
	}
}

template <typename T>
inline T RandomRange(T min, T max) {
    T scale = rand() / (T)RAND_MAX;
    return min + scale * (max - min);
}

std::pair<double, double> AnalyticalSolutionCoeffs(const CartGrid& grid)
{
    auto it = 0;
    double r_inner = 0.1;
    double r_outer = 0.4;
    double bc_inner = 2.0;
    double bc_outer = 6.0;
    for (const auto [key, boundary] : grid.GetImmersedBoundaries())
    {
        if (it == 0)
        {
            r_inner = boundary->GetSize();
            r_outer = boundary->GetSize();
        }
        else
        {
            r_inner = std::min(r_inner, boundary->GetSize());
            r_outer = std::max(r_outer, boundary->GetSize());
        }
        it++;
    }

    for (const auto [key, boundary] : grid.GetImmersedBoundaries())
    {
        if (r_inner == boundary->GetSize()) // dangerous?
        {
            bc_inner = boundary->GetBoundaryPhi();
        }
        if (r_outer == boundary->GetSize())
        {
            bc_outer = boundary->GetBoundaryPhi();
        }
    }

    // phi = A log(r) + B 

    auto a = (bc_inner - bc_outer) / (std::log(r_inner) - std::log(r_outer));
    auto b = bc_inner - a*std::log(r_inner);

    /*std::vector<double> xs(100);
    std::vector<double> ys(100);
    std::generate(xs.begin(), xs.end(), [n = 0]() mutable { return 0.01 * n++; });

    for (size_t i = 0; i < xs.size(); i++)
    {
        ys[i] = a * std::log(xs[i]) + b;
    }

    if (ImGui::Begin("Some plot"))
    {
        ImPlot::BeginPlot("Line Plots");
        ImPlot::SetupAxes("x", "y");
        ImPlot::PlotLine("f(x)", xs.data(), ys.data(), xs.size());
        ImPlot::EndPlot();
        ImGui::End();
    }*/

    return { a, b };
}

Eigen::MatrixXd AnalyticalSolution(const CartGrid& grid)
{
    auto coeffs = AnalyticalSolutionCoeffs(grid);

    Eigen::MatrixXd analytical_solution = Eigen::MatrixXd::Zero(grid.GetPhiMatrix().rows(), grid.GetPhiMatrix().cols());
    //std::cout << "ANALYTICAL\n" << analytical_solution << "\n\n";

    // phi = A log(r) + B 
    for (size_t i = 0; i < grid.GetPhiMatrix().cols(); i++)
    {
        for (size_t j = 0; j < grid.GetPhiMatrix().rows(); j++)
        {
            if (grid.GetCellFlag(i, j) == 1)
            {
                
            }
            else
            {
                Eigen::Vector2d grid_coordinate{ i, j };

                auto world_coordinate_r = grid.GetWorldCoordinate(grid_coordinate) - Eigen::Vector2d{ 0.5, 0.5 };
                auto r = std::sqrt(std::pow(world_coordinate_r.x(), 2) + std::pow(world_coordinate_r.y(), 2));

                analytical_solution(j, i) = coeffs.first * std::log(r) + coeffs.second;
            }
        }
    }

    return analytical_solution;
}

void DataViewer::DataViewerInitialize()
{
    // set up references
    m_data_export.SetSolverRef(m_solver);

    // Setup SDL
    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) != 0)
    {
        printf("Error: %s\n", SDL_GetError());
    }

    // GL 4.6 + GLSL 460
    const char* glsl_version = "#version 460";
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, 0);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 6);

    // From 2.0.18: Enable native IME.
#ifdef SDL_HINT_IME_SHOW_UI
    SDL_SetHint(SDL_HINT_IME_SHOW_UI, "1");
#endif

    // Create window with graphics context
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
    window = SDL_CreateWindow("Data Visualizer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1800, 1000, window_flags); // 1280 x 720
	gl_context = SDL_GL_CreateContext(window);
    SDL_GL_MakeCurrent(window, gl_context);
    gladLoadGL();
    SDL_GL_SetSwapInterval(1); // Enable vsync

    glEnable(GL_DEBUG_OUTPUT);
    glDebugMessageCallback(MessageCallback, 0);

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;           // Enable Docking
    io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;         // Enable Multi-Viewport / Platform Windows

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();

    // When viewports are enabled we tweak WindowRounding/WindowBg so platform windows can look identical to regular ones.
    ImGuiStyle& style = ImGui::GetStyle();
    if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
    {
        style.WindowRounding = 0.0f;
        style.Colors[ImGuiCol_WindowBg].w = 1.0f;
    }

    // Setup Platform/Renderer backends
    ImGui_ImplSDL2_InitForOpenGL(window, gl_context);
    ImGui_ImplOpenGL3_Init(glsl_version);

    // Load Fonts
    // - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
    // - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
    // - If the file cannot be loaded, the function will return NULL. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
    // - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which ImGui_ImplXXXX_NewFrame below will call.
    // - Use '#define IMGUI_ENABLE_FREETYPE' in your imconfig file to use Freetype for higher quality font rendering.
    // - Read 'docs/FONTS.md' for more instructions and details.
    // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
    // - Our Emscripten build process allows embedding fonts to be accessible at runtime from the "fonts/" folder. See Makefile.emscripten for details.
    //io.Fonts->AddFontDefault();
    //io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\segoeui.ttf", 18.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
    //ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, NULL, io.Fonts->GetGlyphRangesJapanese());
    //IM_ASSERT(font != NULL);
}


void DataViewer::RunDataViewer()
{
    ImGuiIO& io = ImGui::GetIO();

	//(void)io;

	while (!quit)
    {
        // Poll and handle events (inputs, window resize, etc.)
        // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
        // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
        // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
        // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            ImGui_ImplSDL2_ProcessEvent(&event);
            if (event.type == SDL_QUIT)
                quit = SDL_TRUE;
            if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_CLOSE && event.window.windowID == SDL_GetWindowID(window))
                quit = SDL_TRUE;

            //User presses a key
            if (event.type == SDL_KEYDOWN)
            {
                camera.ProcessKeyboardCommands(event.key.keysym.sym);
                //Select surfaces based on key press
                switch (event.key.keysym.sym)
                {
                default:
                    break;
                }
            }
        }

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame();
        ImGui::NewFrame();

        // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
        if (show_demo_window)
            ImGui::ShowDemoWindow(&show_demo_window);

        // Main Menu --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        static ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_None;

        // We are using the ImGuiWindowFlags_NoDocking flag to make the parent window not dockable into,
        // because it would be confusing to have two docking targets within each others.
        ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking;
        
        const ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(viewport->WorkPos);
        ImGui::SetNextWindowSize(viewport->WorkSize);
        ImGui::SetNextWindowViewport(viewport->ID);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
        window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
        window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;

        // When using ImGuiDockNodeFlags_PassthruCentralNode, DockSpace() will render our background
        // and handle the pass-thru hole, so we ask Begin() to not render a background.
        if (dockspace_flags & ImGuiDockNodeFlags_PassthruCentralNode)
            window_flags |= ImGuiWindowFlags_NoBackground;

        // Important: note that we proceed even if Begin() returns false (aka window is collapsed).
        // This is because we want to keep our DockSpace() active. If a DockSpace() is inactive,
        // all active windows docked into it will lose their parent and become undocked.
        // We cannot preserve the docking relationship between an active window and an inactive docking, otherwise
        // any change of dockspace/settings would lead to windows being stuck in limbo and never being visible.

        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
        ImGui::Begin("DockSpace Demo", nullptr, window_flags);
        ImGui::PopStyleVar();
		ImGui::PopStyleVar(2);

        // Submit the DockSpace
        ImGuiIO& io = ImGui::GetIO();
        
        ImGuiID dockspace_id = ImGui::GetID("MyDockSpace");
        ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), dockspace_flags);
        

        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("Options"))
            {
                if (ImGui::MenuItem("Save", ""))
                {
	                if (m_solver)
	                {
                        for (const auto& boundary : m_boundaries)
                        {
                            m_data_export.WriteGeometry(boundary->name_id, *dynamic_cast<Circle2D_SDF*>(boundary.get()), boundary->GetSize());
                        }
                        m_data_export.WriteSteadyStateAll();
	                }
                }
                ImGui::Separator();
                ImGui::EndMenu();
            }

            ImGui::EndMenuBar();
        }

        static bool reset_layout = true;
        if (ImGui::DockBuilderGetNode(dockspace_id) == nullptr || reset_layout)
        {
            reset_layout = false;
            ImGui::DockBuilderRemoveNode(dockspace_id);
            ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), dockspace_flags);

            ImGuiID dockspace_id_down = ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Down, 0.3f, nullptr, &dockspace_id);
            ImGuiID dockspace_id_down_left = ImGui::DockBuilderSplitNode(dockspace_id_down, ImGuiDir_Left, 0.4f, nullptr, &dockspace_id_down);
            ImGuiID dockspace_id_left = ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 0.3f, nullptr, &dockspace_id);
            
            ImGui::DockBuilderDockWindow("CENTER", dockspace_id);
            ImGui::DockBuilderDockWindow("LEFT", dockspace_id_left);

            ImGui::DockBuilderDockWindow("DOWN LEFT", dockspace_id_down_left);
            ImGui::DockBuilderDockWindow("DOWN", dockspace_id_down);
            
            
            

            ImGui::DockBuilderFinish(dockspace_id);
        }

        static int selected_simulation_run;
        if (ImGui::Begin("LEFT"))
        {
            if (ImGui::Button("Add Simulation"))
                ImGui::OpenPopup("Create New Simulation");
            ImGui::SameLine();
            if (ImGui::Button("Refine Selected") && !models.empty())
            {
                size_t size = models[selected_mat].m_size;
                size = m_solver->AddGridDoubledSolution(*m_solver->GetSolution(size)); // add new solution and get the new size

                if (size != 0)
                {
                    SolutionModel view_model{};
                    view_model.SetSolution(m_solver->GetSolution(size));
                    view_model.InitData();
                    models.push_back(view_model);
                    selected_mat++;
                }
            }

            // Always center this window when appearing
            ImVec2 center = ImGui::GetMainViewport()->GetCenter();
            ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

            if (ImGui::BeginPopupModal("Create New Simulation", NULL, ImGuiWindowFlags_AlwaysAutoResize))
            {
                ImGui::Text("Create new simulation with the following parameters:\n  size:  grid size (n x n)\n    dt:  time step\n\n");

                const double  dt_min = 0., dt_max = 1.0;
                const unsigned int  size_min = 2, size_max = 1000;
                static unsigned int size_slider = 2.0;
            	static double dt_slider = 0.0;

                ImGui::PushItemWidth(120);
                ImGui::DragScalar("##off", ImGuiDataType_U32, &size_slider, 0.2f, &size_min, &size_max, "size: %u");
                ImGui::SameLine(140);
                ImGui::DragScalar("##off2", ImGuiDataType_Double, &dt_slider, 0.0005f, &dt_min, &dt_max, "dt: %.10f");
                ImGui::PopItemWidth();

                ImGui::NewLine();
            	ImGui::Separator();

                if (ImGui::Button("OK", ImVec2(120, 0)))
                {
                    m_solver->AddSolution(dt_slider, size_slider);

                	SolutionModel view_model{};
                    view_model.SetSolution(m_solver->GetSolution(size_slider));
                    view_model.InitData();
                    models.push_back(view_model);

	                ImGui::CloseCurrentPopup();
                }
                ImGui::SetItemDefaultFocus();
                ImGui::SameLine();
                if (ImGui::Button("Cancel", ImVec2(120, 0)))
                {
	                ImGui::CloseCurrentPopup();
                }
                ImGui::EndPopup();
            }

            // later in your code...
            if (ImGui::BeginCombo("combo", models.size() > 0 ? models[selected_mat].size_id.c_str() : "##nothing")) {
                for (int i = 0; i < models.size(); ++i) {
                    const bool isSelected = (selected_mat == i);
                    if (ImGui::Selectable(models[i].size_id.c_str(), isSelected)) {
                        selected_mat = i;
                    }

                    // Set the initial focus when opening the combo
                    // (scrolling + keyboard navigation focus)
                    if (isSelected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }

            ImGui::Separator();
            ImGui::Text("Boundaries");

            if (ImGui::Button("Add Boundary"))
                ImGui::OpenPopup("Create New Boundary");
            ImGui::SameLine();
            if (ImGui::Button("Clear All") && !models.empty())
            {
                m_boundaries = {};

                // must clean up more..
            }

            // Always center this window when appearing
            center = ImGui::GetMainViewport()->GetCenter();
            ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

            if (ImGui::BeginPopupModal("Create New Boundary", NULL, ImGuiWindowFlags_AlwaysAutoResize))
            {
                ImGui::SeparatorText("Boundary Type");
                std::array<std::string, 1> boundary_types{ "Circle" };
                static int boundary_type_selector = 0;

                if (ImGui::BeginCombo("Type", boundary_types[boundary_type_selector].c_str())) {
                    for (int i = 0; i < boundary_types.size(); ++i) {
                        const bool isSelected = (boundary_type_selector == i);
                        if (ImGui::Selectable(boundary_types[i].c_str(), isSelected)) {
                            boundary_type_selector = i;
                        }

                        // Set the initial focus when opening the combo
                        // (scrolling + keyboard navigation focus)
                        if (isSelected) {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }
                static std::string default_name_str = "Inner Circle";
                static char name_input[128] = "";
                ImGui::InputTextWithHint("##no_text", default_name_str.c_str(), name_input, IM_ARRAYSIZE(name_input));

                std::array<std::string, 2> normal_dir_str{ "Outward", "Inward" };
                static bool normal_dir = 0;
                const double  r_min = 0., r_max = 1.0;
                static double r_slider = 0.0;
                if (boundary_type_selector == 0)
                {
                    ImGui::SeparatorText("Geometry");
                    ImGui::PushItemWidth(140);
                    ImGui::DragScalar("##radius", ImGuiDataType_Double, &r_slider, 0.005f, &r_min, &r_max, "radius: %.4f");

                    std::string text = "Normal Dir (" + normal_dir_str[normal_dir] + ")";
                    if (ImGui::Checkbox(text.c_str(), &normal_dir))
                    {
                    }
                    (normal_dir) ? default_name_str = "Inner Circle" : default_name_str = "Outer Circle";
                    ImGui::PopItemWidth();
                }

                ImGui::SeparatorText("Boundary Condition");

            	std::array<std::string, 2> bc_types{ "Dirichlet", "Neumann" };
                static int bc_type_selector = 0;

                if (ImGui::BeginCombo("BC", bc_types[bc_type_selector].c_str())) {
                    for (int i = 0; i < bc_types.size(); ++i) {
                        const bool isSelected = (bc_type_selector == i);
                        if (ImGui::Selectable(bc_types[i].c_str(), isSelected)) {
                            bc_type_selector = i;
                        }

                        // Set the initial focus when opening the combo
                        // (scrolling + keyboard navigation focus)
                        if (isSelected) {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }


                const double  bc_min = -10., bc_max = 10.0;
                static double bc_slider = 0.0;

                std::string label = "";
                switch (bc_type_selector)
                {
                case 0:
                    label = "phi: %.3f";
                    break;
                case 1:
                    label = "d/dn(phi): %.3f";
                    break;
                default:
                    label = "%.3f";
                }

                ImGui::PushItemWidth(140);
                ImGui::DragScalar("##none", ImGuiDataType_Double, &bc_slider, 0.05f, &bc_min, &bc_max, label.c_str());
                ImGui::PopItemWidth();

                ImGui::NewLine();
                ImGui::Separator();

                if (ImGui::Button("OK", ImVec2(120, 0)))
                {
                    auto geom = std::make_shared<Circle2D_SDF>( Circle2D_SDF{ 0.5, 0.5, r_slider, bc_slider, normal_dir } );

                    m_boundaries.push_back(geom);
                    m_boundaries.back()->name_id = default_name_str;
                    ImGui::CloseCurrentPopup();
                }
                ImGui::SetItemDefaultFocus();
                ImGui::SameLine();
                if (ImGui::Button("Cancel", ImVec2(120, 0)))
                {
                    ImGui::CloseCurrentPopup();
                }
                ImGui::EndPopup();
            }

            static ImGuiTableFlags flags =
                ImGuiTableFlags_Resizable | ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders
        		| ImGuiTableFlags_ScrollY
                | ImGuiTableFlags_SizingFixedFit;

            if (ImGui::BeginTable("Boundaries", 4, flags, {0, 80}))
            {
                
                // Display headers so we can inspect their interaction with borders.
                // (Headers are not the main purpose of this section of the demo, so we are not elaborating on them too much. See other sections for details)
                
                ImGui::TableSetupColumn("Type");
                ImGui::TableSetupColumn("Size");
                ImGui::TableSetupColumn("Boundary Condition");
                ImGui::TableSetupColumn("BC Value");
                ImGui::TableHeadersRow();

                for (int row = 0; row < m_boundaries.size(); row++)
                {
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0);
                    std::string text = "Circle";
                    ImGui::Text(text.c_str());

                    ImGui::TableSetColumnIndex(1);
                    text = std::to_string(m_boundaries[row]->GetSize());
                    ImGui::Text(text.c_str());

                    ImGui::TableSetColumnIndex(2);
                    switch (m_boundaries[row]->GetBoundaryCondition())
                    {
                    case BoundaryCondition::Dirichlet:
                        text = "Dirichlet";
                        break;
                    case BoundaryCondition::Neumann:
                        text = "Neumann";
                        break;
                    default:
                        text = "error";
                    }
                    ImGui::Text(text.c_str());

                    ImGui::TableSetColumnIndex(3);
                    text = std::to_string(m_boundaries[row]->GetBoundaryPhi());
                    ImGui::Text(text.c_str());
                }
                ImGui::EndTable();
            }

            ImGui::SeparatorText("Simulation");

            if (ImGui::Button("Apply Boundaries"))
            {
                for (auto& [size, solution] : m_solver->GetSolutions())
                {
                    for (auto& boundary : m_boundaries)
                    {
                        solution->m_mesh_grid->AddImmersedBoundary(boundary->name_id, boundary);
                    }
                    solution->m_mesh_grid->UpdateGrid();
                }
            }

            if (ImGui::Button("Initialize"))
            {
                for (auto& [size, solution] : m_solver->GetSolutions())
                {
                    solution->m_mesh_grid->InitializeField();
                }
            }

            ImGui::Separator();
            static ImGuiTableFlags table_flags =
                ImGuiTableFlags_Resizable | ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders
                | ImGuiTableFlags_ScrollY
                | ImGuiTableFlags_SizingFixedFit;

            if (ImGui::BeginTable("Simulations", 3, table_flags, { 0, 80 }))
            {

                // Display headers so we can inspect their interaction with borders.
                // (Headers are not the main purpose of this section of the demo, so we are not elaborating on them too much. See other sections for details)

                ImGui::TableSetupColumn("Size");
                ImGui::TableSetupColumn("Iterations");
                ImGui::TableSetupColumn("Time");
                ImGui::TableHeadersRow();

                for (const auto& [size, solution] : m_solver->GetSolutions())
                {
                    //if (!filter.PassFilter(item->Name))
                    //    continue;

                    const bool item_is_selected = table_selection.contains(size);// (selected_simulation_run == size);
                    ImGui::PushID(size);
                    ImGui::TableNextRow(ImGuiTableRowFlags_None, 0.0f);

                    // For the demo purpose we can select among different type of items submitted in the first column
                    ImGui::TableSetColumnIndex(0);
                    std::string label = std::to_string(size) + " x " + std::to_string(size);
                    if (ImGui::Selectable(label.c_str(), item_is_selected, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowItemOverlap, ImVec2(0, 0)))
                    {
                        if (ImGui::GetIO().KeyCtrl)
                        {
                            if (item_is_selected)
                                table_selection.find_erase_unsorted(size);
                            else
                                table_selection.push_back(size);
                        }
                        else
                        {
                            table_selection.clear();
                            table_selection.push_back(size);
                            selected_simulation_run = size;

                            for (int i = 0; i < models.size(); ++i)
                            {
	                            if (models[i].m_size == size)
	                            {
                                    selected_mat = i;
	                            }
                            }
                        }
                    }
                    

                    //ImGui::TableNextRow();
                    //ImGui::TableSetColumnIndex(0);
                    //ImGui::Text("%zu x %zu", size, size);

                    ImGui::TableSetColumnIndex(1);
                    std::string text = std::to_string(solution->m_iteration);
                    ImGui::Text(text.c_str());

                    ImGui::TableSetColumnIndex(2);
                    ImGui::Text(std::to_string(solution->m_time).c_str());

                    ImGui::PopID();
                }
                ImGui::EndTable();
            }


            // Simulation running logic
            const unsigned int  it_min = 0, it_max = 1000;
            static unsigned int it_slider = 100;

            ImGui::PushItemWidth(120);
            if (ImGui::Button("Run Simulation") && selected_simulation_run != 0)
            {
                m_run_simulation = true;
                m_interations_remaining = it_slider;
            }
            ImGui::SameLine(140);
            ImGui::DragScalar("##iterations_slider", ImGuiDataType_U32, &it_slider, 0.2f, &it_min, &it_max, "n: %u");
            ImGui::PopItemWidth();

            if ( m_run_simulation && m_interations_remaining-- > 0 )
            {
                m_solver->GetSolution(selected_simulation_run)->RecursiveUpdateFromThis();
            }
            else
            {
                m_run_simulation = false;
            }

            ImGui::SeparatorText("Richardson Extrapolation");

            if (ImGui::Button("Create Group (Current Selection)"))
            {
                for (const size_t size : table_selection)
                {
                    m_re_group.try_emplace(size, m_solver, m_solver->GetSolution(selected_simulation_run));
                }
            }


            if (ImGui::BeginTable("Solutions", 3, table_flags, { 0, 80 }))
            {

                // Display headers so we can inspect their interaction with borders.
                // (Headers are not the main purpose of this section of the demo, so we are not elaborating on them too much. See other sections for details)

                ImGui::TableSetupColumn("Size");
                ImGui::TableSetupColumn("Iterations");
                ImGui::TableSetupColumn("Time");
                ImGui::TableHeadersRow();

                for (const auto& [size, solution] : m_re_group)
                {
                    //if (!filter.PassFilter(item->Name))
                    //    continue;

                    const bool item_is_selected = table_selection.contains(size);// (selected_simulation_run == size);
                    ImGui::PushID(size);
                    ImGui::TableNextRow(ImGuiTableRowFlags_None, 0.0f);

                    // For the demo purpose we can select among different type of items submitted in the first column
                    ImGui::TableSetColumnIndex(0);
                    std::string label = std::to_string(size) + " x " + std::to_string(size);
                    if (ImGui::Selectable(label.c_str(), item_is_selected, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowItemOverlap, ImVec2(0, 0)))
                    {
                        if (ImGui::GetIO().KeyCtrl)
                        {
                            /*if (item_is_selected)
                                table_selection.find_erase_unsorted(size);
                            else
                                table_selection.push_back(size);*/
                        }
                        else
                        {
                            table_selection.clear();
                            table_selection.push_back(size);
                            selected_simulation_run = size;

                            for (int i = 0; i < models.size(); ++i)
                            {
                                if (models[i].m_size == size)
                                {
                                    selected_mat = i;
                                }
                            }
                        }
                    }


                    //ImGui::TableNextRow();
                    //ImGui::TableSetColumnIndex(0);
                    //ImGui::Text("%zu x %zu", size, size);

                    ImGui::TableSetColumnIndex(1);
                    std::string text = std::to_string(solution.m_coarse_solution->m_time);
                    ImGui::Text(text.c_str());

                    ImGui::TableSetColumnIndex(2);
                    ImGui::Text(std::to_string(solution.m_fine_solution->m_time).c_str());

                    ImGui::PopID();
                }
                ImGui::EndTable();
            }


            static unsigned int re_it_slider = 1;

            ImGui::PushItemWidth(120);
            if (ImGui::Button("Run RE Iteration") && !m_re_group.empty())
            {
                m_re_group.at(selected_simulation_run).Update();
                m_data_export.WriteRichardsonExtrapolationData(m_re_group.at(selected_simulation_run), selected_simulation_run);
            }
            ImGui::SameLine(140);
            ImGui::DragScalar("##re_iterations_slider", ImGuiDataType_U32, &re_it_slider, 0.2f, &it_min, &it_max, "n: %u");
            ImGui::PopItemWidth();


            ImGui::End();
        }

        if (ImGui::Begin("CENTER"))
        {
            {
                ImGui::BeginChild("ChildL", ImVec2(ImGui::GetContentRegionAvail().x * 0.5f, ImGui::GetContentRegionAvail().y), true, NULL);

                static float values1[7][7] = { {0.8f, 2.4f, 2.5f, 3.9f, 0.0f, 4.0f, 0.0f},
                                    {2.4f, 0.0f, 4.0f, 1.0f, 2.7f, 0.0f, 0.0f},
                                    {1.1f, 2.4f, 0.8f, 4.3f, 1.9f, 4.4f, 0.0f},
                                    {0.6f, 0.0f, 0.3f, 0.0f, 3.1f, 0.0f, 0.0f},
                                    {0.7f, 1.7f, 0.6f, 2.6f, 2.2f, 6.2f, 0.0f},
                                    {1.3f, 1.2f, 0.0f, 0.0f, 0.0f, 3.2f, 5.1f},
                                    {0.1f, 2.0f, 0.0f, 1.4f, 0.0f, 1.9f, 6.3f} };

                int size = models[selected_mat].m_size;
                Eigen::MatrixXd analytical = AnalyticalSolution(*m_solver->GetSolution(size)->m_mesh_grid);
                Eigen::MatrixXd error = m_solver->GetSolution(size)->m_mesh_grid->GetPhiMatrixRef() - analytical;
                double* values = error.data();
                //std::cout << "ERROR\n" << error << "\n\n";

                //double* values = m_solver->GetSolution(size)->m_mesh_grid->GetPhiMatrixRef().data();

                static float scale_min = -0.15f;
                static float scale_max = 0.15f;

                static bool auto_scaling = true;
                if (auto_scaling)
                {
                    scale_min = std::max((float)error.minCoeff(), -std::numeric_limits<float>::max());
                    scale_max = std::min((float)error.maxCoeff(), std::numeric_limits<float>::max());
                }

                static const char* xlabels[] = { "C1","C2","C3","C4","C5","C6","C7" };
                static const char* ylabels[] = { "R1","R2","R3","R4","R5","R6","R7" };

                int standard_width = (ImGui::GetContentRegionAvail().x > 225) ? 225 : ImGui::GetContentRegionAvail().x;
                static ImPlotColormap map = ImPlotColormap_Viridis;
                if (ImPlot::ColormapButton(ImPlot::GetColormapName(map), ImVec2(ImGui::GetContentRegionAvail().x, 0), map)) {
                    map = (map + 1) % ImPlot::GetColormapCount();
                    // We bust the color cache of our plots so that item colors will
                    // resample the new colormap in the event that they have already
                    // been created. See documentation in implot.h.
                    ImPlot::BustColorCache("##Heatmap1");
                    ImPlot::BustColorCache("##Heatmap2");
                }

                ImGui::SetNextItemWidth(standard_width);
                ImGui::DragFloatRange2("##min_max", &scale_min, &scale_max, 0.0001f, -10, 10, "Min %.2g", "Max %.2g");


                ImGui::SetNextItemWidth(standard_width);
                ImGui::Checkbox("Colorbar Autoscale", &auto_scaling);
                ImGui::SameLine();
                static ImPlotHeatmapFlags hm_flags = 0;
                ImGui::CheckboxFlags("Column Major", (unsigned int*)&hm_flags, ImPlotHeatmapFlags_ColMajor);

                static ImPlotAxisFlags axes_flags = ImPlotAxisFlags_Lock | ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoTickMarks;

                ImPlot::PushColormap(map);

                double zoom = 0.0;
                static int digits = 0;

                int color_bar_width = 80;
                int map_size = ImGui::GetContentRegionAvail().x - color_bar_width;
                if (ImPlot::BeginPlot("##Heatmap1", ImVec2(map_size, map_size), ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText)) {
                    //ImPlot::SetupAxes(NULL, NULL, axes_flags, axes_flags);
                    //ImPlot::SetupAxisTicks(ImAxis_X1, 0 + 1.0 / 14.0, 1 - 1.0 / 14.0, 7, xlabels);
                    //ImPlot::SetupAxisTicks(ImAxis_Y1, 1 - 1.0 / 14.0, 0 + 1.0 / 14.0, 7, ylabels);
                    //ImPlot::PlotHeatmap("heat", values1[0], 7, 7, scale_min, scale_max, "%g", ImPlotPoint(0, 0), ImPlotPoint(1, 1), hm_flags);

                    ImPlot::SetupAxes(NULL, NULL, ImPlotAxisFlags_NoDecorations, ImPlotAxisFlags_NoDecorations);
                    ImPlot::SetupAxesLimits(0, 1, 0, 1);

                    ImPlotContext& gp = *ImPlot::GetCurrentContext();
                    zoom = gp.CurrentPlot->Axes->Range.Size();

                    digits = std::clamp(static_cast<int>(0.372 * std::pow(zoom / (31.3248 * std::pow(size, -1.00715)), -1.554)), 0, 8);
                    double treshold = 31.3248 * std::pow(size, -1.00715);
                    std::string format = (treshold < zoom) ? "" : "%." + std::to_string(digits) + "f";

                    ImPlot::PlotHeatmap("##heat1", values, size, size, scale_min, scale_max, format.c_str());

                	ImPlot::EndPlot();
                }
                ImGui::SameLine();
                ImPlot::ColormapScale("##HeatScale", scale_min, scale_max, ImVec2(color_bar_width, map_size), "%.0e");
                
                ImPlot::PopColormap();

                ImGui::Text("%g", zoom);

                ImGui::EndChild();
            }

	        if (!models.empty())
	        {
                ImGui::SameLine();
                ImGui::BeginChild("TextureRender");
                // Get the size of the child (i.e. the whole draw size of the windows).
                ImVec2 wsize = ImGui::GetWindowSize();
                ImVec2 wpos = ImGui::GetWindowPos();
                ImVec2 mouse_pos = ImGui::GetMousePos();

                glm::vec2 rel_mouse_pos { std::clamp((mouse_pos.x - wpos.x) / wsize.x, 0.0f, 1.0f), std::clamp((mouse_pos.y - wpos.y) / wsize.y, 0.0f, 1.0f) };
                static glm::vec2 rel_mouse_pos_old = rel_mouse_pos;
                static glm::vec2 delta_old{ 0.0 };
                static glm::vec2 delta{ 0.0 };

                if (io.MouseDownDuration[0] < 0.01f) // just pressed
                {
	                // reset mouse pos
                    rel_mouse_pos_old = rel_mouse_pos;
                }
                delta_old = delta;
                delta = rel_mouse_pos - rel_mouse_pos_old;

                auto acc = delta - delta_old;
				acc.x = -acc.x;

                if (io.MouseDownDuration[0] > 0.01f && ImGui::IsWindowFocused())
                {
                    camera.MouseTranslate(acc);
                }
                if (ImGui::IsWindowHovered())
                {
                    camera.ScrollZoom(io.MouseWheel);
                }

                // Because I use the texture from OpenGL, I need to invert the V from the UV.
                ImTextureID id = (ImTextureID)models[selected_mat].GetTextureBuffer();
                //ImTextureID id = (ImTextureID)test_meshes[selected_mat].out_texture_buffer;
                ImGui::Image(id, wsize, ImVec2(0, 1), ImVec2(1, 0));
                ImGui::EndChild();
	        }
            
            ImGui::End();
        }

        if (ImGui::Begin("DOWN"))
        {
            const double  scaling_min = 0.5, scaling_max = 1.5;
            static double a_d_scaling = 1.0;

            static unsigned int wlsq_idx = 0;
            unsigned int size = 0;
            double a_d = 0.0;

            ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x/2 * 0.75);
            ImGui::DragScalar("##a_d_scaling", ImGuiDataType_Double, &a_d_scaling, 0.05f, &scaling_min, &scaling_max, "scaling factor: %g");
            ImGui::SameLine();
            if(ImGui::Button("Update") && selected_simulation_run != 0)
            {
                m_solver->GetSolution(selected_simulation_run)->m_mesh_grid->m_weight_scaling = a_d_scaling;
                m_solver->GetSolution(selected_simulation_run)->m_mesh_grid->WLSQUpdateGeometry();
            }

            ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x/2);
            ImGui::DragScalar("##wlsq_gp_select", ImGuiDataType_U32, &wlsq_idx, 0.2f, 0, &size-1, "idx: %u");

            if (ImPlot::BeginPlot("Weighting Function", ImVec2{ ImGui::GetContentRegionAvail().x/2, ImGui::GetContentRegionAvail().y })) {
                ImPlot::SetupAxes("distance", "weight");
                //ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
                
                if (selected_simulation_run != 0 )
                {
                    auto solution = m_solver->GetSolution(selected_simulation_run);
                    if ( !solution->m_mesh_grid->GetWLSQdata().empty() )
                    {
                        auto& wlsq = solution->m_mesh_grid->GetWLSQdata();
                        auto keys = solution->m_mesh_grid->GetGhostPoints();
                        size = keys.size();

                        wlsq_idx = std::clamp(wlsq_idx, (unsigned int)0, size-1);

                        auto& selected_wlsq = wlsq.at(keys.at(wlsq_idx));
 
                        a_d = selected_wlsq.m_a_d;

                        ImPlot::PlotScatter("##active_points", selected_wlsq.dist.data(), selected_wlsq.weight.data(), selected_wlsq.dist.size());
                    }
                }

                auto res = 1000;
                std::vector<double> xs(res);
                std::vector<double> ys(res);
                std::generate(xs.begin(), xs.end(), [n = 0, &ys, a_d, res]() mutable
                {
                	auto x = 50.0/res * n++;
                    ys[n] = std::exp(-std::pow( x, 2.0 ) / a_d);

	                return x;
                });
            	
            	ImPlot::PlotLine("##no_labels", xs.data(), ys.data(), xs.size());

                ImPlot::EndPlot();
            }

            
            ImGui::End();
        }

        if (ImGui::Begin("DOWN LEFT"))
        {
            ImPlot::ShowDemoWindow();

            if (ImPlot::BeginPlot("Convergence", ImVec2{ ImGui::GetContentRegionAvail().x, ImGui::GetContentRegionAvail().y})) {
                ImPlot::SetupAxes("time", "L2");
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                //ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
                for (const auto& [key, solution] : m_solver->GetSolutions())
                {
                    auto dt = solution->m_dt;
                    std::vector<double> t(solution->m_l2_norms.size());
                    std::generate(t.begin(), t.end(), [n = 0, &dt]() mutable { return dt * n++; });

                    ImPlot::PlotLine(std::to_string(key).c_str(), t.data(), solution->m_l2_norms.data(), solution->m_l2_norms.size());
                }
                
                ImPlot::EndPlot();
            }

            ImGui::End();
        }


        ImGui::End();
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // Rendering
        ImGui::Render();
        glViewport(0, 0, (int)io.DisplaySize.x, (int)io.DisplaySize.y);
        glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);

        // OpenGL stuff

        models[selected_mat].SetMVP(camera.GetMVP());
        models[selected_mat].DrawSolutionModelToTexture();

        // ImGui
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        // Update and Render additional Platform Windows
        // (Platform functions may change the current OpenGL context, so we save/restore it to make it easier to paste this code elsewhere.
        //  For this specific demo app we could also call SDL_GL_MakeCurrent(window, gl_context) directly)
        if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
        {
            SDL_Window* backup_current_window = SDL_GL_GetCurrentWindow();
            SDL_GLContext backup_current_context = SDL_GL_GetCurrentContext();
            ImGui::UpdatePlatformWindows();
            ImGui::RenderPlatformWindowsDefault();
            SDL_GL_MakeCurrent(backup_current_window, backup_current_context);
        }


        SDL_GL_SwapWindow(window);
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();

    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();
}