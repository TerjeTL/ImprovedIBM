
#pragma once

#include "ibm_application/SDLGraphics.h"
#include <iomanip>



SDLGraphics::SDLGraphics(std::shared_ptr<CartGrid> mesh_grid) : m_mesh_grid(mesh_grid)
{
    auto grid_extent = m_mesh_grid->GetMeshSize();
    grid_width = grid_extent.first;
    grid_height = grid_extent.second;

    window_width = (grid_width * grid_cell_size) + 2;
    window_height = (grid_height * grid_cell_size) + 2;

    grid_cursor = {
        ((float)grid_width - 1.f) / 2.f * (float)grid_cell_size,
        ((float)grid_height - 1.f) / 2.f * (float)grid_cell_size,
        (float)grid_cell_size,
        (float)grid_cell_size
    };

    grid_position = {
        (grid_width - 1) / 2,
        (grid_height - 1) / 2
    };

    grid_cursor_ghost = { 0, 0, grid_cell_size, grid_cell_size };
}

void SDLGraphics::SDLGraphicsInitialize()
{
    SDL_SetMainReady();
    if (SDL_Init(SDL_INIT_EVERYTHING) < 0)
    {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Initialize SDL: %s",
            SDL_GetError());
    }
    TTF_Init();
    
    if (SDL_CreateWindowAndRenderer(window_width, window_height, SDL_WINDOW_SHOWN | SDL_RENDERER_PRESENTVSYNC, &window, &renderer) < 0)
    {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION,
            "Create window and renderer: %s", SDL_GetError());
    }
    font = TTF_OpenFont("C:/Development/C++/ImprovedIBM/fonts/Roboto/Roboto-Regular.ttf", 64);

    SDL_SetWindowTitle(window, "SDL IBM");
}

ScreenSpacePos SDLGraphics::GetScreenSpacePos(GridPos grid_location)
{
    int even_grid_offset_x = 0;
    int even_grid_offset_y = 0;
    even_grid_offset_x = grid_cell_size / 2;
    even_grid_offset_y = grid_cell_size / 2;
    
    int x_screenspace = (grid_location.x - grid_position.x  + grid_width / 2)*grid_cell_size + grid_cell_size/2 - zoom_level * grid_width / 2 - even_grid_offset_x;
    int y_screenspace = (grid_location.y - grid_position.y  + grid_height / 2)*grid_cell_size + grid_cell_size/2 - zoom_level * grid_height / 2 - even_grid_offset_y;
    return { x_screenspace, y_screenspace };
}

ScreenSpacePosF SDLGraphics::GetScreenSpacePosF(GridPosF grid_location)
{
    double even_grid_offset_x = 0.0;
    double even_grid_offset_y = 0.0;
    even_grid_offset_x = 0.5 * static_cast<double>(grid_cell_size);
    even_grid_offset_y = 0.5 * static_cast<double>(grid_cell_size);
    
    double x_screenspace = (grid_location.x - static_cast<double>(grid_position.x) + static_cast<double>(grid_width) / 2.0) * static_cast<double>(grid_cell_size)
        + static_cast<double>(grid_cell_size) / 2.0
        - static_cast<double>(zoom_level) * static_cast<double>(grid_width) / 2.0
        - even_grid_offset_x;
    double y_screenspace = (grid_location.y - static_cast<double>(grid_position.y) + static_cast<double>(grid_height) / 2.0) * static_cast<double>(grid_cell_size)
        + static_cast<double>(grid_cell_size) / 2.0
        - static_cast<double>(zoom_level) * static_cast<double>(grid_height) / 2.0
        - even_grid_offset_y;
    return { x_screenspace, y_screenspace };
}

GridPosF SDLGraphics::GetGridPosF(double x_pos_norm, double y_pos_norm)
{
    double x_gridspace = (static_cast<double>(grid_width)-1) * x_pos_norm;
    double y_gridspace = (static_cast<double>(grid_height)-1) * y_pos_norm;

    return { x_gridspace, y_gridspace };
}

void SDLGraphics::RenderCircle(SDL_Renderer* renderer, int32_t centreX, int32_t centreY, int32_t radius, SDL_Color color)
{
    SDL_Color original;
    SDL_GetRenderDrawColor(renderer, &original.r, &original.g, &original.b, &original.a);
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);

    const int32_t diameter = (radius * 2);

    int32_t x = (radius - 1);
    int32_t y = 0;
    int32_t tx = 1;
    int32_t ty = 1;
    int32_t error = (tx - diameter);

    while (x >= y)
    {
        //  Each of the following renders an octant of the circle
        SDL_RenderDrawPoint(renderer, centreX + x, centreY - y);
        SDL_RenderDrawPoint(renderer, centreX + x, centreY + y);
        SDL_RenderDrawPoint(renderer, centreX - x, centreY - y);
        SDL_RenderDrawPoint(renderer, centreX - x, centreY + y);
        SDL_RenderDrawPoint(renderer, centreX + y, centreY - x);
        SDL_RenderDrawPoint(renderer, centreX + y, centreY + x);
        SDL_RenderDrawPoint(renderer, centreX - y, centreY - x);
        SDL_RenderDrawPoint(renderer, centreX - y, centreY + x);

        if (error <= 0)
        {
            ++y;
            error += ty;
            ty += 2;
        }

        if (error > 0)
        {
            --x;
            tx += 2;
            error += (tx - diameter);
        }
    }

    SDL_SetRenderDrawColor(renderer, original.r, original.g, original.b, original.a);
}

void SDLGraphics::RenderFillCircle(SDL_Renderer* renderer, int32_t centreX, int32_t centreY, int32_t radius, SDL_Color color)
{
    SDL_Color original;
    SDL_GetRenderDrawColor(renderer, &original.r, &original.g, &original.b, &original.a);
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);

    const int32_t diameter = (radius * 2);

    int32_t x = (radius - 1);
    int32_t y = 0;
    int32_t tx = 1;
    int32_t ty = 1;
    int32_t error = (tx - diameter);

    while (x >= y)
    {
        //  Each of the following renders an octant of the circle
        SDL_RenderDrawLine(renderer, centreX + x, centreY + y, centreX - x, centreY + y);
        SDL_RenderDrawLine(renderer, centreX + y, centreY + x, centreX - y, centreY + x);
        SDL_RenderDrawLine(renderer, centreX + x, centreY - y, centreX - x, centreY - y);
        SDL_RenderDrawLine(renderer, centreX + y, centreY - x, centreX - y, centreY - x);

        if (error <= 0)
        {
            ++y;
            error += ty;
            ty += 2;
        }

        if (error > 0)
        {
            --x;
            tx += 2;
            error += (tx - diameter);
        }
    }

    SDL_SetRenderDrawColor(renderer, original.r, original.g, original.b, original.a);
}

void SDLGraphics::RenderText(std::string text, GridPosF position, SDL_Color color, float width, float height)
{
    m_message_surface = TTF_RenderText_Solid(font, text.c_str(), color);

    // Convert to texture
    m_message = SDL_CreateTextureFromSurface(renderer, m_message_surface);

    ScreenSpacePosF number_pos = GetScreenSpacePosF(position);

    SDL_FRect message_rect;
    message_rect.x = number_pos.x;
    message_rect.y = number_pos.y;
    message_rect.w = width;
    message_rect.h = height;

    SDL_RenderCopyF(renderer, m_message, NULL, &message_rect);

    SDL_FreeSurface(m_message_surface);
    SDL_DestroyTexture(m_message);
}

void SDLGraphics::RenderGridText()
{
    if (grid_cell_size / 2 < 30 || !m_draw_node_values)
    {
        return;
    }

    for (int i = 0; i < grid_width; i += 1)
    {
        for (int j = 0; j < grid_height; j += 1)
        {
            SDL_Color text_color = { 255, 255, 255 };

            // Get phi value as string
            std::ostringstream value_ostream;
            value_ostream << std::showpoint << std::setprecision(4) << m_mesh_grid->GetPhi(i,j);
            std::string value_str = value_ostream.str();

            RenderText(value_str, GridPosF{ (double)i, (double)j - 0.2 }, text_color, 40.0, 15.0);

            if (m_mesh_grid->GetCellFlag(i, j) == 2)
            {
                value_ostream.str("");
                value_ostream.clear();

                value_ostream << std::showpoint << std::setprecision(4) << m_mesh_grid->GetImagePointPhi(i,j);
                value_str = value_ostream.str();

                Eigen::Vector2d image_point_pos = m_mesh_grid->GetGridCoordinate(m_mesh_grid->GetImagePoint(i, j));

                GridPosF image_point_text_pos = { image_point_pos(0), image_point_pos(1) - 0.15 };

                RenderText(value_str, image_point_text_pos, text_color, 40.0, 15.0);
            }

        }
    }
}

void SDLGraphics::RunSDLGraphics()
{
    while (!quit) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
            case SDL_KEYDOWN:
                switch (event.key.keysym.sym) {
                case SDLK_w:
                case SDLK_UP:
                    grid_position.y -= 1;
                    std::cout << "P = (" << grid_position.x << ", " << grid_position.y << ")\n";
                    //std::cout << immersed_boundary.SignedDistanceFunction((double)grid_position.x, (double)grid_position.y) << "\n";
                    break;
                case SDLK_s:
                case SDLK_DOWN:
                    grid_position.y += 1;
                    //std::cout << immersed_boundary.SignedDistanceFunction((double)grid_position.x, (double)grid_position.y) << "\n";
                    std::cout << "P = (" << grid_position.x << ", " << grid_position.y << ")\n";
                    break;
                case SDLK_a:
                case SDLK_LEFT:
                    grid_position.x -= 1;
                    //std::cout << immersed_boundary.SignedDistanceFunction((double)grid_position.x, (double)grid_position.y) << "\n";
                    std::cout << "P = (" << grid_position.x << ", " << grid_position.y << ")\n";
                    break;
                case SDLK_d:
                case SDLK_RIGHT:
                    grid_position.x += 1;
                    //std::cout << immersed_boundary.SignedDistanceFunction((double)grid_position.x, (double)grid_position.y) << "\n";
                    std::cout << "P = (" << grid_position.x << ", " << grid_position.y << ")\n";
                    break;
                case SDLK_z:
                    if (zoom_gain <= 0)
                    {
                        zoom_gain = 1;
                    }
                    grid_cell_size += zoom_gain;
                    zoom_level += zoom_gain;
                    zoom_gain += zoom_gain;
                    zoom_gain = std::min(zoom_gain, 32);
                    std::cout << "zoom " << zoom_level << "\n";
                    break;
                case SDLK_x:
                    if (zoom_gain >= 0)
                    {
                        zoom_gain = -1;
                    }
                  
                    grid_cell_size += zoom_gain;
                    
                    zoom_level += zoom_gain;
                    zoom_gain += zoom_gain;
                    zoom_gain = std::max(zoom_gain, -16);
                    
                    if (grid_cell_size < 2)
                    {
                        auto diff = -2 - grid_cell_size;
                        grid_cell_size += diff;
                        zoom_level += diff;
                        zoom_gain = -1;
                        break;
                    }
                    
                    std::cout << "zoom " << zoom_level << "\n";
                    break;
                case SDLK_v:
                    m_draw_node_values = !m_draw_node_values;
                    break;
                }
                break;
            case SDL_MOUSEBUTTONDOWN:
                //grid_position.x = (event.motion.x / grid_cell_size) * grid_cell_size;
                //grid_position.y = (event.motion.y / grid_cell_size) * grid_cell_size;
                break;
            case SDL_MOUSEMOTION:
                //grid_cursor_ghost.x = (event.motion.x / grid_cell_size) * grid_cell_size;
                //grid_cursor_ghost.y = (event.motion.y / grid_cell_size) * grid_cell_size;

                if (!mouse_active)
                    mouse_active = SDL_TRUE;
                break;
            case SDL_WINDOWEVENT:
                if (event.window.event == SDL_WINDOWEVENT_ENTER && !mouse_hover)
                    mouse_hover = SDL_TRUE;
                else if (event.window.event == SDL_WINDOWEVENT_LEAVE && mouse_hover)
                    mouse_hover = SDL_FALSE;
                break;
            case SDL_QUIT:
                quit = SDL_TRUE;
                break;
            }
        }

        float center_cursor_size = (float)grid_cell_size / 4.f;
        ScreenSpacePosF center_cursor_pos{ (double)window_width * 0.5, (double)window_height * 0.5 };
        center_cursor_pos.x -= 0.5*center_cursor_size;
        center_cursor_pos.y -= 0.5*center_cursor_size;

        grid_cursor = {
        (float)center_cursor_pos.x,
        (float)center_cursor_pos.y,
        center_cursor_size,
        center_cursor_size
        };

        // Draw grid background.
        SDL_SetRenderDrawColor(renderer, grid_background.r, grid_background.g,
            grid_background.b, grid_background.a);
        SDL_RenderClear(renderer);

        // Draw grid lines.
        SDL_SetRenderDrawColor(renderer, grid_line_color.r, grid_line_color.g,
            grid_line_color.b, grid_line_color.a);

        ScreenSpacePos west_north_lim = GetScreenSpacePos(GridPos{ 0, 0 });
        ScreenSpacePos east_south_lim = GetScreenSpacePos(GridPos{ grid_width, grid_height });

        int west_lim = west_north_lim.x;
        int east_lim = east_south_lim.x;

        int north_lim = west_north_lim.y;
        int south_lim = east_south_lim.y;

        for (int x = 0; x < grid_width; x += 1)
        {
            ScreenSpacePos screenspace_pos = GetScreenSpacePos(GridPos{ x, 0 });
            SDL_RenderDrawLine(renderer, screenspace_pos.x, south_lim - grid_cell_size, screenspace_pos.x, north_lim);
            for (int y = 0; y < grid_height; y += 1)
            {
                screenspace_pos = GetScreenSpacePos(GridPos{ x, y });
                SDL_RenderDrawLine(renderer, west_lim, screenspace_pos.y, east_lim - grid_cell_size, screenspace_pos.y);
            }
        }

        // Draw our grid flags
        if (m_mesh_grid)
        {
            for (int i = 0; i < grid_width; i += 1)
            {
                for (int j = 0; j < grid_height; j += 1)
                {
                    ScreenSpacePos screenspace_pos = GetScreenSpacePos(GridPos{ i, j });

                    switch (m_mesh_grid->GetCellFlag(i,j))
                    {
                    case 0:
                        RenderFillCircle(renderer, screenspace_pos.x, screenspace_pos.y, grid_cell_size / 20);
                        break;
                    case 1:
                        RenderFillCircle(renderer, screenspace_pos.x, screenspace_pos.y, grid_cell_size / 20, SDL_Color{ 255, 0, 0, 255 });
                        break;
                    case 2:
                    {
                        RenderFillCircle(renderer, screenspace_pos.x, screenspace_pos.y, grid_cell_size / 20, SDL_Color{ 0, 255, 0, 255 });
  
                        auto image_point = m_mesh_grid->GetImagePoint(i, j);
                        ScreenSpacePosF image_point_screenspace = GetScreenSpacePosF(GetGridPosF(image_point[0], image_point[1]));
                        RenderCircle(renderer, (float)image_point_screenspace.x, (float)image_point_screenspace.y, grid_cell_size / 20, SDL_Color{ 0, 0, 255, 255 });
                        break;
                    }
                    default:
                        break;
                    }
                }
            }
        }

        GridPosF debug_point_grid_from = GetGridPosF(0.0, 0.0);
        ScreenSpacePosF debug_point_scr_from = GetScreenSpacePosF(debug_point_grid_from);

        GridPosF debug_point_grid_to = GetGridPosF(0.5, 0.5);
        ScreenSpacePosF debug_point_scr_to = GetScreenSpacePosF(debug_point_grid_to);
        
        SDL_RenderDrawLineF(renderer, debug_point_scr_from.x, debug_point_scr_from.y, debug_point_scr_to.x, debug_point_scr_to.y);
        RenderFillCircle(renderer, debug_point_scr_to.x, debug_point_scr_to.y, grid_cell_size / 20);

        // Draw grid ghost cursor.
        if (mouse_active && mouse_hover) {
            SDL_SetRenderDrawColor(renderer, grid_cursor_ghost_color.r,
                grid_cursor_ghost_color.g,
                grid_cursor_ghost_color.b,
                grid_cursor_ghost_color.a);
            SDL_RenderFillRect(renderer, &grid_cursor_ghost);
        }

        // Draw grid cursor.
        SDL_SetRenderDrawColor(renderer, grid_cursor_color.r,
            grid_cursor_color.g, grid_cursor_color.b,
            grid_cursor_color.a);
        SDL_RenderFillRectF(renderer, &grid_cursor);


        //ScreenSpacePos inner_circle_center = GetScreenSpacePos(GridPos{ 10, 10 });
        //RenderFillCircle(renderer, inner_circle_center.x, inner_circle_center.y, grid_cell_size*4);

        for (auto const& [hash, sdf] : m_mesh_grid->GetImmersedBoundaries())
        {
            sdf->RenderSDF(renderer, *this);
        }

        RenderGridText();

        SDL_RenderPresent(renderer);
    }

    TTF_CloseFont(font);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    TTF_Quit();
    SDL_Quit();
}