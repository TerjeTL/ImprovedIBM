#pragma once

#include "data_viewer/DataViewer.h"
#include "data_viewer/Shader.h"
#include "data_viewer/Camera.h"
#include "data_viewer/SolutionModel.h"

#include "glad/glad.h"
#include "imgui.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_opengl3.h"
#include <stdio.h>
#include <SDL.h>
#include <SDL_opengl.h>
#include <vector>
#include "glm.hpp"
#include "imgui_internal.h"

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
    window = SDL_CreateWindow("Data Visualizer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1280, 720, window_flags);
	gl_context = SDL_GL_CreateContext(window);
    SDL_GL_MakeCurrent(window, gl_context);
    gladLoadGL();
    SDL_GL_SetSwapInterval(1); // Enable vsync

    glEnable(GL_DEBUG_OUTPUT);
    glDebugMessageCallback(MessageCallback, 0);

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
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
        static bool opt_fullscreen = true;
        static ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_None;

        // We are using the ImGuiWindowFlags_NoDocking flag to make the parent window not dockable into,
        // because it would be confusing to have two docking targets within each others.
        ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking;
        if (opt_fullscreen)
        {
            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImGui::SetNextWindowPos(viewport->WorkPos);
            ImGui::SetNextWindowSize(viewport->WorkSize);
            ImGui::SetNextWindowViewport(viewport->ID);
            ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
            ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
            window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
            window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;
        }
        else
        {
            dockspace_flags &= ~ImGuiDockNodeFlags_PassthruCentralNode;
        }

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
        if (opt_fullscreen)
            ImGui::PopStyleVar(2);

        // Submit the DockSpace
        ImGuiIO& io = ImGui::GetIO();
        
        ImGuiID dockspace_id = ImGui::GetID("MyDockSpace");
        ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), dockspace_flags);
        

        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("Options"))
            {
                // Disabling fullscreen would allow the window to be moved to the front of other windows,
                // which we can't undo at the moment without finer window depth/z control.
                ImGui::MenuItem("Fullscreen", NULL, &opt_fullscreen);
                ImGui::Separator();

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

            ImGuiID dockspace_id_down = ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Down, 0.25f, nullptr, &dockspace_id);
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

                    const bool item_is_selected = (selected_simulation_run == size);
                    ImGui::PushID(size);
                    ImGui::TableNextRow(ImGuiTableRowFlags_None, 0.0f);

                    // For the demo purpose we can select among different type of items submitted in the first column
                    ImGui::TableSetColumnIndex(0);
                    std::string label = std::to_string(size) + " x " + std::to_string(size);
                    if (ImGui::Selectable(label.c_str(), item_is_selected, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowItemOverlap, ImVec2(0, 0)))
                    {
                        if (ImGui::GetIO().KeyCtrl)
                        {
                            //if (item_is_selected)
                            //    selection.find_erase_unsorted(size);
                            //else
                            //    selection.push_back(size);
                        }
                        else
                        {
                            selected_simulation_run = size;
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
                }
                ImGui::EndTable();
            }


            // Simulation running logic
            const unsigned int  it_min = 0, it_max = 1000;
            static unsigned int it_slider = 800;

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
            
            ImGui::End();
        }

        if (ImGui::Begin("CENTER"))
        {
            {
                ImGui::BeginChild("ChildL", ImVec2(ImGui::GetContentRegionAvail().x * 0.25f, ImGui::GetContentRegionAvail().y), true, NULL);
                for (int i = 0; i < 10; i++)
                    ImGui::Text("%04d: scrollable region", i);
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
            
            ImGui::End();
        }

        if (ImGui::Begin("DOWN LEFT"))
        {

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

        for (auto& model : models)
        {
            model.SetMVP(camera.GetMVP());
            model.DrawSolutionModelToTexture();
        }

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
    ImGui::DestroyContext();

    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();
}