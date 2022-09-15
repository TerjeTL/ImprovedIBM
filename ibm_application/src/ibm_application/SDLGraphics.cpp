
#pragma once

#include "ibm_application/SDLGraphics.h"



SDLGraphics::SDLGraphics()
{
}

void SDLGraphics::SDLGraphicsInitialize()
{
    SDL_SetMainReady();
    if (SDL_Init(SDL_INIT_EVERYTHING) < 0)
    {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Initialize SDL: %s",
            SDL_GetError());
    }

    if (SDL_CreateWindowAndRenderer(window_width, window_height, SDL_WINDOW_SHOWN | SDL_RENDERER_PRESENTVSYNC, &window, &renderer) < 0)
    {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION,
            "Create window and renderer: %s", SDL_GetError());
    }

    SDL_SetWindowTitle(window, "SDL IBM");
}

ScreenSpacePos SDLGraphics::GetScreenSpacePos(GridPos grid_location)
{
    int even_grid_offset_x = 0;
    int even_grid_offset_y = 0;
    if (grid_width % 2 == 0)
    {
        even_grid_offset_x = grid_cell_size / 2;
    }
    if (grid_height % 2 == 0)
    {
        even_grid_offset_y = grid_cell_size / 2;
    }
    int x_screenspace = (grid_location.x - grid_position.x  + grid_width / 2)*grid_cell_size + grid_cell_size/2 - zoom_level * grid_width / 2 - even_grid_offset_x;
    int y_screenspace = (grid_location.y - grid_position.y  + grid_height / 2)*grid_cell_size + grid_cell_size/2 - zoom_level * grid_height / 2 - even_grid_offset_y;
    return { x_screenspace, y_screenspace };
}

ScreenSpacePosF SDLGraphics::GetScreenSpacePosF(GridPosF grid_location)
{
    double even_grid_offset_x = 0.0;
    double even_grid_offset_y = 0.0;
    if (grid_width % 2 == 0)
    {
        even_grid_offset_x = 0.5 * static_cast<double>(grid_cell_size);
    }
    if (grid_height % 2 == 0)
    {
        even_grid_offset_y = 0.5 * static_cast<double>(grid_cell_size);
    }
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
    double x_gridspace = static_cast<double>(grid_width * grid_cell_size) * x_pos_norm;
    double y_gridspace = static_cast<double>(grid_height * grid_cell_size) * y_pos_norm;

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

        for (int x = 0; x < grid_width; x += 1) 
        {
            for (int y = 0; y < grid_height; y += 1)
            {
                ScreenSpacePos screenspace_pos = GetScreenSpacePos(GridPos{ x, y });
                
                if (immersed_boundary.SignedDistanceFunction(x, y) < 0)
                {
                    if (immersed_boundary.SignedDistanceFunction(x + 1, y) >= 0
                        || immersed_boundary.SignedDistanceFunction(x + -1, y) >= 0
                        || immersed_boundary.SignedDistanceFunction(x, y+1) >= 0
                        || immersed_boundary.SignedDistanceFunction(x, y-1) >= 0)
                    {
                        RenderFillCircle(renderer, screenspace_pos.x, screenspace_pos.y, grid_cell_size / 20, SDL_Color{ 0, 255, 0, 255 });
                    }
                    else
                    {
                        RenderFillCircle(renderer, screenspace_pos.x, screenspace_pos.y, grid_cell_size / 20, SDL_Color{ 255, 0, 0, 255 });
                    }
                }
                else
                {
                    RenderFillCircle(renderer, screenspace_pos.x, screenspace_pos.y, grid_cell_size / 20);
                }
            }
        }

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

        immersed_boundary.RenderSDF(renderer, *this);

        //ScreenSpacePos outer_circle_center = GetScreenSpacePos(GridPos{ grid_width/2, grid_height/2 });
        //RenderCircle(renderer, outer_circle_center.x, outer_circle_center.y, grid_cell_size * grid_height/6);

        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}