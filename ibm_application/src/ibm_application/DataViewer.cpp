#pragma once

#include "ibm_application/DataViewer.h"
#include "data_viewer/Shader.h"
#include "data_viewer/Camera.h"

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

void GLAPIENTRY
MessageCallback(GLenum source,
    GLenum type,
    GLuint id,
    GLenum severity,
    GLsizei length,
    const GLchar* message,
    const void* userParam)
{
    fprintf(stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
        (type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""),
        type, severity, message);
}


template<typename T>
std::vector<T> EigenMatrixToArray(const Eigen::MatrixX<T>& eigen_mat)
{
    std::vector<T> mat_vector(eigen_mat.rows()*eigen_mat.cols());

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajMat;
    RowMajMat::Map(mat_vector.data(), eigen_mat.rows(), eigen_mat.cols()) = eigen_mat;

    return mat_vector;
}

class MatrixModel {
    unsigned int texture1, texture2, outTexture;
    unsigned int VBO, VAO, EBO;
    const float vertices[20] = {
        //pos                //texcoords
         0.5f,  0.5f, 0.0f,  1.0f, 1.0f, //tr
         0.5f, -0.5f, 0.0f,  1.0f, 0.0f, //br
        -0.5f, -0.5f, 0.0f,  0.0f, 0.0f, //bl
        -0.5f,  0.5f, 0.0f,  0.0f, 1.0f  //tl
    };
    const unsigned int indices[6] = {
        0, 1, 3, // first triangle
        1, 2, 3  // second triangle
    };

    glm::mat4 MVP;
    glm::vec3 color;

    Eigen::Matrix4f test_tex = (Eigen::Matrix4f() <<
        0.0f, 0.05f, 0.1f, 0.15f,
        0.2f, 0.25f, 0.3f, 0.35f,
        0.4f, 0.45f, 0.5f, 0.55f,
        0.6f, 0.65f, 0.7f, 1.0f ).finished();

    Shader shaders{ "D:/Development/C++/ibm_application/ibm_application/src/shaders/MatrixViewer.vert",
	"D:/Development/C++/ibm_application/ibm_application/src/shaders/MatrixViewer.frag" };
public:
    unsigned int framebuffer, textureColorbuffer;
    std::vector<int> grid_flag_mat;
    uint32_t size = 0;

    MatrixModel(std::vector<int> grid_flags) : grid_flag_mat(grid_flags)
	{
        color = glm::vec3(1, 0, 0);
        MVP = glm::mat4(1.0f);

        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);

        glBindVertexArray(VAO);

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), &vertices[0], GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), &indices[0], GL_STATIC_DRAW);

        // position attribute
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        // texture coord attribute
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(1);

        glBindVertexArray(0);

        // load and create a texture 
		// -------------------------
        // texture 1 (float)
        // ---------
        glGenTextures(1, &texture1);
        glBindTexture(GL_TEXTURE_2D, texture1);
        // set the texture wrapping parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);	// set texture wrapping to GL_REPEAT (default wrapping method)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        float borderColor[] = { 1.0f, 1.0f, 0.0f, 1.0f };
        glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);
        // set texture filtering parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    	//float* data = test_tex.data();
        auto vec_data = EigenMatrixToArray<float>(test_tex);
    	float* data = vec_data.data();
        if (data)
        {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 4, 4, 0, GL_RED, GL_FLOAT, data);
            glGenerateMipmap(GL_TEXTURE_2D);
        }
        else
        {
            std::cout << "Failed to load texture" << std::endl;
        }

        // load and create a texture 
        // -------------------------
        // texture 2 (int)
        // ---------
        glGenTextures(1, &texture2);
        glBindTexture(GL_TEXTURE_2D, texture2);
        // set the texture wrapping parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);	// set texture wrapping to GL_REPEAT (default wrapping method)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);
        // set texture filtering parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

        int* data_2 = grid_flag_mat.data();
    	size = std::round(std::sqrt(grid_flag_mat.size()));
        if (data_2)
        {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_R16I, size, size, 0, GL_RED_INTEGER, GL_INT, data_2);
            glGenerateMipmap(GL_TEXTURE_2D);
        }
        else
        {
            std::cout << "Failed to load texture" << std::endl;
        }

        // tell opengl for each sampler to which texture unit it belongs to (only has to be done once)
        // -------------------------------------------------------------------------------------------
        shaders.use(); // don't forget to activate/use the shader before setting uniforms!
        // either set it manually like so:
        // or set it via the texture class
        shaders.setFloat("texture1", 0);
        shaders.setInt("texture2", 1);

        glGenFramebuffers(1, &framebuffer);
        glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
        glGenTextures(1, &textureColorbuffer);
        glBindTexture(GL_TEXTURE_2D, textureColorbuffer);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 800, 800, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureColorbuffer, 0);

        // now that we actually created the framebuffer and added all attachments we want to check if it is actually complete now
        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
            std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!" << std::endl;
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    void SetMVP(glm::mat4 mvp) {
        MVP = mvp;
    }

    int setColor(glm::vec3 color_in) {
        color = color_in;
        return 1;
    }

    int draw() {
        shaders.use();
        
    	// bind textures on corresponding texture units
        glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        glViewport(0, 0, 800, 800);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture1);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, texture2);

        shaders.SetMat4fv("mvp", MVP);
        shaders.SetVec3fv("color", color);
        shaders.SetUInt("size", size);

        glBindVertexArray(VAO);

        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        glBindVertexArray(0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        return 1;
    }

    /*~MatrixModel() {

        glDeleteVertexArrays(1, &VAO);
        glDeleteBuffers(1, &VBO);
        glDeleteBuffers(1, &EBO);
    }*/
};

void DataViewer::DataViewerInitialize()
{
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
    Camera camera;

    std::vector<MatrixModel> models;

    models.emplace_back(EigenMatrixToArray<int>(grids[0]->GetGridFlags()));
    models.emplace_back(EigenMatrixToArray<int>(grids[1]->GetGridFlags()));

    ImGuiIO& io = ImGui::GetIO(); (void)io;
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
        static bool opt_padding = false;
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
                ImGui::MenuItem("Padding", NULL, &opt_padding);
                ImGui::Separator();

                if (ImGui::MenuItem("Flag: NoSplit", "", (dockspace_flags & ImGuiDockNodeFlags_NoSplit) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_NoSplit; }
                if (ImGui::MenuItem("Flag: NoResize", "", (dockspace_flags & ImGuiDockNodeFlags_NoResize) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_NoResize; }
                if (ImGui::MenuItem("Flag: NoDockingInCentralNode", "", (dockspace_flags & ImGuiDockNodeFlags_NoDockingInCentralNode) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_NoDockingInCentralNode; }
                if (ImGui::MenuItem("Flag: AutoHideTabBar", "", (dockspace_flags & ImGuiDockNodeFlags_AutoHideTabBar) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_AutoHideTabBar; }
                if (ImGui::MenuItem("Flag: PassthruCentralNode", "", (dockspace_flags & ImGuiDockNodeFlags_PassthruCentralNode) != 0, opt_fullscreen)) { dockspace_flags ^= ImGuiDockNodeFlags_PassthruCentralNode; }
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

        if (ImGui::Begin("LEFT"))
        {
            std::vector<std::string> mats;
            for (const auto& mat : grids)
            {
                auto size = mat->GetMeshSize();
                std::string str = std::to_string(size.first) + "x" + std::to_string(size.second);
                mats.push_back(str);

            }

            // later in your code...
            if (ImGui::BeginCombo("combo", mats[selected_mat].c_str())) {
                for (int i = 0; i < mats.size(); ++i) {
                    const bool isSelected = (selected_mat == i);
                    if (ImGui::Selectable(mats[i].c_str(), isSelected)) {
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

            ImGui::End();
        }

        if (ImGui::Begin("CENTER"))
        {
            ImGui::BeginChild("GameRender");
            // Get the size of the child (i.e. the whole draw size of the windows).
            ImVec2 wsize = ImGui::GetWindowSize();
            // Because I use the texture from OpenGL, I need to invert the V from the UV.
            ImTextureID id = (ImTextureID)models[selected_mat].textureColorbuffer;
            ImGui::Image(id, wsize, ImVec2(0, 1), ImVec2(1, 0));
            ImGui::EndChild();
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
        //line_1.draw();

        for (auto& model : models)
        {
            model.SetMVP(camera.GetMVP());
            model.draw();
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