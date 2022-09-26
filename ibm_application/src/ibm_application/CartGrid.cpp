
#include "ibm_application/CartGrid.h"

#include <iostream>

CartGrid::CartGrid(size_t nn, double phi_0) : grid_flags(Eigen::MatrixXi::Zero(nn, nn)), phi_matrix(Eigen::MatrixXd::Zero(nn, nn)),
phi_image_point_matrix(Eigen::MatrixXd::Zero(nn, nn)), boundary_phi(Eigen::MatrixXd::Zero(nn, nn))
{
    h = length_scales(0) / static_cast<double>(nn - 1);
}

void CartGrid::UpdateGrid()
{
    for (auto const& [hash, sdf] : immersed_boundaries)
    {
        // For each boundary (given that the boundary data exists) we update the grid flags
        if (sdf)
        {
            for (size_t i = 0UL; i < grid_flags.rows(); ++i) {
                for (size_t j = 0UL; j < grid_flags.cols(); ++j) {

                    auto x = static_cast<double>(i) * h;
                    auto y = static_cast<double>(j) * h;

                    // Check if i,j should be a ghost, active or inactive grid point and assign correct flag
                    if (sdf->SignedDistanceFunction(x, y) < 0.0) {
                        if (sdf->SignedDistanceFunction(x + h, y) >= 0.0
                            || sdf->SignedDistanceFunction(x - h, y) >= 0.0
                            || sdf->SignedDistanceFunction(x, y + h) >= 0.0
                            || sdf->SignedDistanceFunction(x, y - h) >= 0.0) {
                            grid_flags(i, j) = 2;

                            // Calculate normal and assign respective image-point (they are not cleared!)
                            auto normal = sdf->GetNormal(x, y);
                            Eigen::Vector2d node{ x, y };

                            Eigen::Vector2d image_point = node + normal * std::abs(sdf->SignedDistanceFunction(x, y)) * 2.0;
                            image_points[{ i, j }] = image_point;

                            phi_image_point_matrix(i,j) = BilinearInterpolation(i, j);
                            boundary_phi(i, j) = sdf->GetBoundaryPhi();
                        }
                        else {
                            grid_flags(i, j) = 1;
                        }
                    }
                }
            }
        }
    }
}

double CartGrid::BilinearInterpolation(size_t i, size_t j)
{
    // Interpolation in grid-space
    Eigen::Vector2d grid_coordinate = GetGridCoordinate(image_points.at({i, j}));

    std::array<Eigen::Vector2i, 4> node_indices;
    std::array<Eigen::Vector2d, 4> node_locations;
    Eigen::Vector4d node_phi;

    node_indices[0] = grid_coordinate.cast<int>();
    node_indices[1] = node_indices[0] + Eigen::Vector2i{ 0,1 };
    node_indices[2] = node_indices[0] + Eigen::Vector2i{ 1,1 };
    node_indices[3] = node_indices[0] + Eigen::Vector2i{ 1,0 };

    for (size_t p = 0; p < node_indices.size(); p++)
    {
        if (grid_flags(node_indices[p](0), node_indices[p](1)) == 2)
        {
            auto image_pt_loc = GetGridCoordinate(GetImagePoint(node_indices[p](0), node_indices[p](1)));

            node_locations[p] = node_indices[p].cast<double>() + (image_pt_loc - node_indices[p].cast<double>())/2.0;
            node_phi(p) = boundary_phi(node_indices[p](0), node_indices[p](1));
        }
        else
        {
            node_locations[p] = node_indices[p].cast<double>();
            node_phi(p) = phi_matrix(node_indices[p](0), node_indices[p](1));
        }
    }

    Eigen::Matrix4d vandermonde{
        {1, node_locations[0](0), node_locations[0](1), node_locations[0](0) * node_locations[0](1)},
        {1, node_locations[1](0), node_locations[1](1), node_locations[1](0) * node_locations[1](1)},
        {1, node_locations[2](0), node_locations[2](1), node_locations[2](0) * node_locations[2](1)},
        {1, node_locations[3](0), node_locations[3](1), node_locations[3](0) * node_locations[3](1)}
    };

    //Eigen::Vector4d coefficients = vandermonde.inverse() * node_phi;
    Eigen::Vector4d coefficients = vandermonde.colPivHouseholderQr().solve(node_phi);
    // solve as linear system x = A\b matlab
    // thomas algorithm/tdma

    double phi = coefficients(0)
        + coefficients(1) * grid_coordinate(0)
        + coefficients(2) * grid_coordinate(1)
        + coefficients(3) * grid_coordinate(0) * grid_coordinate(1);

    bilinear_interp_selection[{i, j}] = node_locations;

    return phi;
}
