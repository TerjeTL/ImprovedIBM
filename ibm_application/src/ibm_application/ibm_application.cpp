// gnuplot_application.cpp : Defines the entry point for the application.
//
#pragma once

#include <stdexcept>
#include <filesystem>
#include <iostream>
#include <fstream>

#include "ibm_application/ibm_application.h"
#include "ibm_application/SDLGraphics.h"
#include "ibm_application/CartGrid.h"
#include "ibm_application/Solver.h"
#include <Eigen/Core>

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

void writeToCSVfile(std::string name, Eigen::MatrixXd matrix)
{
    std::filesystem::path file_path = std::filesystem::current_path() / name;
    std::ofstream file(file_path.string());
    file << matrix.format(CSVFormat);
    file.close();
}

int main()
{
    // Debugging grid
    std::shared_ptr<CartGrid> grid_debug = std::make_shared<CartGrid>(22);
    Solver test_solver{ 0.00001, std::make_unique<FTCS_Scheme>(grid_debug), grid_debug };

    // Prepare an SDLGraphics instance
    SDLGraphics sdl_program(grid_debug);
   
    bool exit_menu = false;
    while (!exit_menu)
    {
        std::cout << "--- Immersed Boundary Method Program ---\n"
            << "Run options:\n"
            << "1. Grid Visualization (SDL)\n"
            << "2. Run one step\n"
            << "3. Save Data\n"
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
            grid_debug->AddImmersedBoundary("Inner Cylinder", std::make_shared<Circle2D_SDF>(Circle2D_SDF{ 0.5, 0.5, 0.15, 100.0 }));
            grid_debug->AddImmersedBoundary("Outer Cylinder", std::make_shared<Circle2D_SDF>(Circle2D_SDF{ 0.5, 0.5, 0.44, 200.0, true }));
            grid_debug->UpdateGrid();

            test_solver.PerformStep(100);
            break;
        }
        case 3:
        {
            std::string filename = "";
            std::cout << "Save data to file:";
            std::cin >> filename;
            std::cout << std::endl;

            writeToCSVfile(filename, grid_debug->GetPhiMatrix());
            break;
        }
        default:
            exit_menu = true;
            break;
        }
    }

    return 0;
}