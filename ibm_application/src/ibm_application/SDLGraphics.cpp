
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
    int x_screenspace = (grid_location.x - grid_position.x  + grid_width / 2)*grid_cell_size + grid_cell_size/2 - zoom_level * grid_width / 2;
    int y_screenspace = (grid_location.y - grid_position.y  + grid_width / 2)*grid_cell_size + grid_cell_size/2 - zoom_level * grid_height / 2;
    return { x_screenspace, y_screenspace };
}

void SDLGraphics::DrawCircle(SDL_Renderer* renderer, int32_t centreX, int32_t centreY, int32_t radius)
{
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
                    break;
                case SDLK_s:
                case SDLK_DOWN:
                    grid_position.y += 1;
                    std::cout << "P = (" << grid_position.x << ", " << grid_position.y << ")\n";
                    break;
                case SDLK_a:
                case SDLK_LEFT:
                    grid_position.x -= 1;
                    std::cout << "P = (" << grid_position.x << ", " << grid_position.y << ")\n";
                    break;
                case SDLK_d:
                case SDLK_RIGHT:
                    grid_position.x += 1;
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

        grid_cursor = {
        ((float)grid_width) / 2.f * (float)grid_cell_size - 0.125f * (float)grid_cell_size - (float)zoom_level * (float)grid_width / 2.f,
        ((float)grid_height) / 2.f * (float)grid_cell_size - 0.125f * (float)grid_cell_size - (float)zoom_level * (float)grid_height / 2.f,
        (float)grid_cell_size/4.f,
        (float)grid_cell_size/4.f
        };

        // Draw grid background.
        SDL_SetRenderDrawColor(renderer, grid_background.r, grid_background.g,
            grid_background.b, grid_background.a);
        SDL_RenderClear(renderer);

        // Draw grid lines.
        SDL_SetRenderDrawColor(renderer, grid_line_color.r, grid_line_color.g,
            grid_line_color.b, grid_line_color.a);

        int west_lim = (grid_width / 2 * grid_cell_size - grid_position.x * grid_cell_size + grid_cell_size / 2 - zoom_level * grid_width/2);
        int east_lim = grid_width * grid_cell_size + west_lim;

        int north_lim = (grid_height / 2 * grid_cell_size - grid_position.y*grid_cell_size + grid_cell_size / 2 - zoom_level * grid_height/2);
        int south_lim = grid_height * grid_cell_size + north_lim;
        //std::cout << left_lim << "\n";

        for (int x = west_lim; x < east_lim; x += grid_cell_size) {
            SDL_RenderDrawLine(renderer, x, south_lim-grid_cell_size, x, north_lim);
        }

        for (int y = north_lim; y < south_lim; y += grid_cell_size) {
            SDL_RenderDrawLine(renderer, west_lim, y, east_lim-grid_cell_size, y);
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


        ScreenSpacePos inner_circle_center = GetScreenSpacePos(GridPos{ 10, 10 });
        DrawCircle(renderer, inner_circle_center.x, inner_circle_center.y, grid_cell_size*4);

        ScreenSpacePos outer_circle_center = GetScreenSpacePos(GridPos{ grid_width/2, grid_height/2 });
        DrawCircle(renderer, outer_circle_center.x, outer_circle_center.y, grid_cell_size * grid_height/2);

        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}