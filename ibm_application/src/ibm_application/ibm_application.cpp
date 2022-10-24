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
#include "ibm_application/RichardsonMethod.h"
#include "ibm_application/DataExporter.h"
#include <Eigen/Core>
#include <highfive/H5Easy.hpp>


const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

void writeToCSVfile(std::string name, Eigen::MatrixXd matrix)
{
    std::filesystem::path file_path = std::filesystem::absolute(name);
    std::ofstream file(file_path.string());
    file << matrix.format(CSVFormat);
    file.close();
}

int main(int argc, char* argv[])
{
    std::shared_ptr<CartGrid> coarse_grid = std::make_shared<CartGrid>(42);
    std::shared_ptr<CartGrid> fine_grid = std::make_shared<CartGrid>(83);
    std::shared_ptr<CartGrid> fine_2 = std::make_shared<CartGrid>(165);
    std::shared_ptr<CartGrid> fine_3 = std::make_shared<CartGrid>(329);
    //std::shared_ptr<CartGrid> fine_4 = std::make_shared<CartGrid>(657);

    std::shared_ptr<RichardsonMethod> richardson_extrapolation = std::make_shared<RichardsonMethod>();
    richardson_extrapolation->AddMeshGrid(0, coarse_grid);
    richardson_extrapolation->AddMeshGrid(1, fine_grid);

    std::shared_ptr<Solver> test_solver = std::make_shared<Solver>(0.0001);
    test_solver->AddSolution(0, std::make_unique<FTCS_Scheme>(coarse_grid), coarse_grid);
    test_solver->AddSolution(1, std::make_unique<FTCS_Scheme>(fine_grid), fine_grid);
    test_solver->AddSolution(2, std::make_unique<FTCS_Scheme>(fine_2), fine_2);
    test_solver->AddSolution(3, std::make_unique<FTCS_Scheme>(fine_3), fine_3);
    //test_solver->AddSolution(4, std::make_unique<FTCS_Scheme>(fine_4), fine_4);
    //test_solver->SetRichardsonMethod(richardson_extrapolation);

    std::shared_ptr<DataExporter> data_export = std::make_shared<DataExporter>(std::filesystem::current_path().parent_path().parent_path() / "scripts/export_data.h5");

    test_solver->SetDataExporter(data_export);

    // Prepare an SDLGraphics instance
    SDLGraphics sdl_program(coarse_grid);
   
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
            std::shared_ptr<Circle2D_SDF> inner_circle = std::make_shared<Circle2D_SDF>(Circle2D_SDF{ 0.5, 0.5, 0.15, 1.0 });
            //inner_circle->SetBoundaryCondition(BoundaryCondition::Neumann);
            
            std::shared_ptr<Circle2D_SDF> outer_circle = std::make_shared<Circle2D_SDF>(Circle2D_SDF{ 0.5, 0.5, 0.45, 2.0, true });

            data_export->WriteGeometry("inner", *inner_circle, 0.15);
            data_export->WriteGeometry("outer", *outer_circle, 0.45);

            coarse_grid->AddImmersedBoundary("Inner Cylinder", inner_circle);
            coarse_grid->AddImmersedBoundary("Outer Cylinder", outer_circle);
            coarse_grid->UpdateGrid();

            fine_grid->AddImmersedBoundary("Inner Cylinder", inner_circle);
            fine_grid->AddImmersedBoundary("Outer Cylinder", outer_circle);
            fine_grid->UpdateGrid();

            fine_2->AddImmersedBoundary("Inner Cylinder", inner_circle);
            fine_2->AddImmersedBoundary("Outer Cylinder", outer_circle);
            fine_2->UpdateGrid();

            fine_3->AddImmersedBoundary("Inner Cylinder", inner_circle);
            fine_3->AddImmersedBoundary("Outer Cylinder", outer_circle);
            fine_3->UpdateGrid();

            //fine_4->AddImmersedBoundary("Inner Cylinder", inner_circle);
            //fine_4->AddImmersedBoundary("Outer Cylinder", outer_circle);
            //fine_4->UpdateGrid();

            test_solver->PerformStep(-1);
            break;
        }
        case 3:
        {
            std::string filename = "phi_matrix.csv";
            
            /*std::cout << "Save data to filename:";
            std::cin >> filename;
            std::cout << std::endl;*/
            
            Eigen::MatrixXd result = coarse_grid->GetPhiMatrix();

            for (size_t i = 0; i < result.cols(); i++)
            {
                for (size_t j = 0; j < result.rows(); j++)
                {
                    if (coarse_grid->GetCellFlag(i, j) != 0)
                    {
                        result(i, j) = 0;
                    }
                }
            }

            if (argc > 1)
            {
                filename = "\\" + filename;
                filename = argv[argc - 1] + filename;
                writeToCSVfile(filename, result);
            }
            break;
        }
        default:
            exit_menu = true;
            break;
        }
    }

    return 0;
}