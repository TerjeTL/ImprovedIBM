﻿// gnuplot_application.cpp : Defines the entry point for the application.
//
#pragma once

#include <stdexcept>
#include <iostream>

#include "ibm_application/ibm_application.h"
#include "ibm_application/SDLGraphics.h"
#include "ibm_application/CartGrid.h"

int main()
{
    // Debugging grid
    std::shared_ptr<CartGrid> grid_debug = std::make_shared<CartGrid>(31);

    // Prepare an SDLGraphics instance
    SDLGraphics sdl_program(grid_debug);
   
    bool exit_menu = false;
    while (!exit_menu)
    {
        std::cout << "--- Immersed Boundary Method Program ---\n"
            << "Run options:\n"
            << "1. Grid Visualization (SDL)\n"
            << "2. Generate CartGrid\n"
            << std::endl;

        std::string input;
        getline(std::cin, input);
        
        int arg = 0;
        try {
            arg = std::stoi(input);
        }
        catch (std::invalid_argument& e) 
        {
            std::cerr << e.what() << std::endl;
            return -1;
        }
     
        switch (arg)
        {
        case 1:
        {
            // Launch our SDLGraphics program
            sdl_program.SDLGraphicsInitialize();
            sdl_program.RunSDLGraphics();
            break;
        }
        case 2:
        {
            grid_debug->AddImmersedBoundary("Inner Cylinder", std::make_shared<Circle2D_SDF>(Circle2D_SDF{ 0.5, 0.5, 0.25 }));
            grid_debug->UpdateGrid();
            break;
        }
        default:
            exit_menu = true;
            break;
        }
    }

    return 0;
}