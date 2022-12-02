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

int EquivalentIterations(int refinement_level, int base_iterations)
{
    int iteration_level = std::pow(2, refinement_level * 2);
    return base_iterations * iteration_level;
}

int main(int argc, char* argv[])
{
    std::shared_ptr<CartGrid> grid_0 = std::make_shared<CartGrid>(12);
    std::shared_ptr<CartGrid> grid_1 = std::make_shared<CartGrid>(23);
    std::shared_ptr<CartGrid> grid_2 = std::make_shared<CartGrid>(45);
    std::shared_ptr<CartGrid> grid_3 = std::make_shared<CartGrid>(89);
    std::shared_ptr<CartGrid> grid_4 = std::make_shared<CartGrid>(177);
    std::shared_ptr<CartGrid> grid_5 = std::make_shared<CartGrid>(353);
    //std::shared_ptr<CartGrid> grid_6 = std::make_shared<CartGrid>(705);
    //std::shared_ptr<CartGrid> grid_7 = std::make_shared<CartGrid>(1409);

    int base_level_iterations = 200;
    std::shared_ptr<Solver> test_solver = std::make_shared<Solver>(0.002, base_level_iterations);
    test_solver->AddSolution(0, std::make_unique<FTCS_Scheme>(grid_0), grid_0, base_level_iterations);
    test_solver->AddSolution(1, std::make_unique<FTCS_Scheme>(grid_1), grid_1, EquivalentIterations(1, base_level_iterations) );
    test_solver->AddSolution(2, std::make_unique<FTCS_Scheme>(grid_2), grid_2, EquivalentIterations(2, base_level_iterations) );
    test_solver->AddSolution(3, std::make_unique<FTCS_Scheme>(grid_3), grid_3, EquivalentIterations(3, base_level_iterations) );
    test_solver->AddSolution(4, std::make_unique<FTCS_Scheme>(grid_4), grid_4, EquivalentIterations(4, base_level_iterations) );
    test_solver->AddSolution(5, std::make_unique<FTCS_Scheme>(grid_5), grid_5, EquivalentIterations(5, base_level_iterations));
    //test_solver->AddSolution(6, std::make_unique<FTCS_Scheme>(grid_6), grid_6, EquivalentIterations(6, base_level_iterations));
    //test_solver->AddSolution(7, std::make_unique<FTCS_Scheme>(grid_7), grid_7, EquivalentIterations(7, base_level_iterations));
    //test_solver->AddSolution(4, std::make_unique<FTCS_Scheme>(fine_4), fine_4);
    //test_solver->SetRichardsonMethod(richardson_extrapolation);

    std::shared_ptr<DataExporter> data_export = std::make_shared<DataExporter>(std::filesystem::current_path().parent_path().parent_path() / "scripts/export_data.h5", DataExporter::LoggingConfig::Steady);

    test_solver->SetDataExporter(data_export);

    // Prepare an SDLGraphics instance
    SDLGraphics sdl_program(grid_1);
   
    bool exit_menu = false;
    while (!exit_menu)
    {
        std::cout << "--- Immersed Boundary Method Program ---\n"
            << "Run options:\n"
            << "1. Grid Visualization (SDL)\n"
            << "2. Run Default Configuration\n"
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

            grid_0->AddImmersedBoundary("Inner Cylinder", inner_circle);
            grid_0->AddImmersedBoundary("Outer Cylinder", outer_circle);
            grid_0->UpdateGrid();
            grid_0->InitializeField();

            grid_1->AddImmersedBoundary("Inner Cylinder", inner_circle);
            grid_1->AddImmersedBoundary("Outer Cylinder", outer_circle);
            grid_1->UpdateGrid();
            grid_1->InitializeField();

            grid_2->AddImmersedBoundary("Inner Cylinder", inner_circle);
            grid_2->AddImmersedBoundary("Outer Cylinder", outer_circle);
            grid_2->UpdateGrid();
            grid_2->InitializeField();

            grid_3->AddImmersedBoundary("Inner Cylinder", inner_circle);
            grid_3->AddImmersedBoundary("Outer Cylinder", outer_circle);
            grid_3->UpdateGrid();
            grid_3->InitializeField();

            grid_4->AddImmersedBoundary("Inner Cylinder", inner_circle);
            grid_4->AddImmersedBoundary("Outer Cylinder", outer_circle);
            grid_4->UpdateGrid();
            grid_4->InitializeField();

            grid_5->AddImmersedBoundary("Inner Cylinder", inner_circle);
            grid_5->AddImmersedBoundary("Outer Cylinder", outer_circle);
            grid_5->UpdateGrid();
            grid_5->InitializeField();

            //grid_6->AddImmersedBoundary("Inner Cylinder", inner_circle);
            //grid_6->AddImmersedBoundary("Outer Cylinder", outer_circle);
            //grid_6->UpdateGrid();
            //grid_6->InitializeField();

            //grid_7->AddImmersedBoundary("Inner Cylinder", inner_circle);
            //grid_7->AddImmersedBoundary("Outer Cylinder", outer_circle);
            //grid_7->UpdateGrid();
            //grid_7->InitializeField();

            //fine_4->AddImmersedBoundary("Inner Cylinder", inner_circle);
            //fine_4->AddImmersedBoundary("Outer Cylinder", outer_circle);
            //fine_4->UpdateGrid();

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
            
            Eigen::MatrixXd result = grid_0->GetPhiMatrix();

            for (size_t i = 0; i < result.cols(); i++)
            {
                for (size_t j = 0; j < result.rows(); j++)
                {
                    if (grid_0->GetCellFlag(i, j) != 0)
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