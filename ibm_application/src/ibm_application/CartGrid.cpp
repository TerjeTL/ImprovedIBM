
// This file defines the cartesian grid

#include "ibm_application/CartGrid.h"

#include <iostream>
#include <assert.h>

CartGrid::CartGrid(size_t nn, double phi_0) : grid_flags(Eigen::MatrixXi::Zero(nn, nn)), phi_matrix(Eigen::MatrixXd::Zero(nn, nn)),
phi_image_point_matrix(Eigen::MatrixXd::Zero(nn, nn)), boundary_phi(Eigen::MatrixXd::Zero(nn, nn)), ghost_point_parent_sdf(Eigen::MatrixX<size_t>::Zero(nn, nn))
{
    h = length_scales(0) / static_cast<double>(nn - 1);
}

// Function to update grid flags and compute beta coefficents
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
                            
                            // Set location descriptors
                            grid_flags(i, j) = 2;
                            ghost_point_parent_sdf(i, j) = hash;

                            // Calculate normal and assign respective image-point (they are not cleared!)
                            auto normal = sdf->GetNormal(x, y);
                            Eigen::Vector2d node{ x, y };

                            Eigen::Vector2d image_point = node + normal * std::abs(sdf->SignedDistanceFunction(x, y)) * m_ip_stencil_length_factor;
                            image_points[{ i, j }] = image_point;

                            // Regenerate coefficients for bilinear interpolation
                            UpdateInterpolationCoeffs(i, j);

                            phi_image_point_matrix(i,j) = BilinearInterpolation(i, j);
                            boundary_phi(i, j) = sdf->GetBoundaryPhi();
                        }
                        else
                        {
                            // Set cell as inactive (solid)
                            grid_flags(i, j) = 1;
                        }
                    }
                }
            }
        }
    }
}


// Apply linear distribution to the initialization
void CartGrid::InitializeField()
{
    assert(immersed_boundaries.size() == 2 && "Can only do custom field initialization with exactly two immersed boundaries");

    for (size_t i = 0UL; i < grid_flags.rows(); ++i) {
        for (size_t j = 0UL; j < grid_flags.cols(); ++j) {

            // Skip inactive cells 
            if (grid_flags(i, j) != 0)
            {
                continue;
            }

            double phi_1 = 0.0;
            double phi_2 = 0.0;

            double dist_1 = 0.0;
            double dist_2 = 0.0;

            auto x = static_cast<double>(i) * h;
            auto y = static_cast<double>(j) * h;

            int k = 1;
            for (auto const& [hash, sdf] : immersed_boundaries)
            {
                double dist = sdf->SignedDistanceFunction(x, y);
                double phi = sdf->GetBoundaryPhi();

                if (k == 1)
                {
                    dist_1 = dist;
                    phi_1 = phi;
                }
                else
                {
                    dist_2 = dist;
                    phi_2 = phi;
                }

                k++;
            }

            // Linearly interpolate phi between boundaries
            phi_matrix(i,j) = dist_1 / std::abs(dist_1 + dist_2) * (phi_2 - phi_1) + phi_1;
        }
    }
}

// Update beta coefficients
void CartGrid::UpdateInterpolationCoeffs(size_t i, size_t j)
{
    // Interpolation in grid-space
    Eigen::Vector2d world_coordinate = image_points.at({ i, j });
    Eigen::Vector2d grid_coordinate = GetGridCoordinate(image_points.at({ i, j }));
    size_t parent_sdf = ghost_point_parent_sdf(i, j);
    BoundaryCondition boundary = immersed_boundaries.at(parent_sdf)->GetBoundaryCondition();

    std::array<Eigen::Vector2i, 4> node_indices;
    std::array<Eigen::Vector2d, 4> node_locations;
    Eigen::Vector4d node_phi;

    node_indices[0] = grid_coordinate.cast<int>();
    node_indices[1] = node_indices[0] + Eigen::Vector2i{ 0,1 };
    node_indices[2] = node_indices[0] + Eigen::Vector2i{ 1,1 };
    node_indices[3] = node_indices[0] + Eigen::Vector2i{ 1,0 };

    Eigen::Matrix4d vandermonde{
        {1, node_locations[0](0), node_locations[0](1), node_locations[0](0) * node_locations[0](1)},
        {1, node_locations[1](0), node_locations[1](1), node_locations[1](0) * node_locations[1](1)},
        {1, node_locations[2](0), node_locations[2](1), node_locations[2](0) * node_locations[2](1)},
        {1, node_locations[3](0), node_locations[3](1), node_locations[3](0) * node_locations[3](1)}
    };

    for (size_t p = 0; p < node_indices.size(); p++)
    {
        Eigen::Vector4d vandermonde_row;

        if (grid_flags(node_indices[p](0), node_indices[p](1)) == 2)
        {
            auto image_pt_loc = GetGridCoordinate(GetImagePoint(node_indices[p](0), node_indices[p](1)));

            auto world_normal = immersed_boundaries.at(parent_sdf)->GetNormal(world_coordinate.x(), world_coordinate.y());
            world_normal.normalize();
            auto grid_normal = GetGridCoordinate(world_normal);

            node_locations[p] = node_indices[p].cast<double>() + (image_pt_loc - node_indices[p].cast<double>()) / 2.0;

            if (boundary == BoundaryCondition::Dirichlet)
            {
                node_phi(p) = boundary_phi(node_indices[p](0), node_indices[p](1));
                vandermonde_row = { 1, node_locations[p](0), node_locations[p](1), node_locations[p](0) * node_locations[p](1) };
            }
            else // Neumann
            {
                node_phi(p) = boundary_phi(node_indices[p](0), node_indices[p](1));
                vandermonde_row = { 0, grid_normal.x(), grid_normal.y(), node_locations[p](0) * grid_normal.y() + node_locations[p](1) * grid_normal.x() };
            }
        }
        else
        {
            node_locations[p] = node_indices[p].cast<double>();
            node_phi(p) = phi_matrix(node_indices[p](0), node_indices[p](1));
            vandermonde_row = { 1, node_locations[p](0), node_locations[p](1), node_locations[p](0) * node_locations[p](1) };
        }

        vandermonde.row(p) = vandermonde_row;
    }

    Eigen::Vector4d ip_coeff_vec = { 1, grid_coordinate.x(), grid_coordinate.y(), grid_coordinate.x() * grid_coordinate.y() };
    beta_coeffs[{i, j}] = vandermonde.inverse().transpose() * ip_coeff_vec;

    bilinear_interp_selection[{i, j}] = node_locations;
}

double CartGrid::BilinearInterpolation(size_t i, size_t j)
{
    // Interpolation in grid-space
    // Eigen::Vector2d world_coordinate = image_points.at({ i, j });
    Eigen::Vector2d grid_coordinate = GetGridCoordinate(image_points.at({i, j}));
    size_t parent_sdf = ghost_point_parent_sdf(i, j);
    BoundaryCondition boundary = immersed_boundaries.at(parent_sdf)->GetBoundaryCondition();

    std::array<Eigen::Vector2i, 4> node_indices;
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
            if (boundary == BoundaryCondition::Dirichlet)
            {
                node_phi(p) = boundary_phi(node_indices[p](0), node_indices[p](1));
            }
            else
            {
                node_phi(p) = boundary_phi(node_indices[p](0), node_indices[p](1));
            }
        }
        else
        {
            node_phi(p) = phi_matrix(node_indices[p](0), node_indices[p](1));
        }
    }

    Eigen::Vector4d coefficients_new = beta_coeffs.at({ i,j });
    double phi = coefficients_new.dot(node_phi);

    return phi;
}

// Implementing WLSQ method

void CartGrid::WeightedLeastSquaresMethod(size_t i, size_t j)
{
    // get the selection of nodes
    
    // start by getting center location of node selection square (god gave us no solution for the circle problem)
    auto image_point = GetGridCoordinate(GetImagePoint(i, j));
    Eigen::Vector2d ghost_point = { static_cast<double>(i), static_cast<double>(j) };

    auto boundary_intercept = ghost_point + (image_point - ghost_point) / 2.0;

    // Find the dimensions needed to grab the required number of points
    int c_x = std::round(boundary_intercept.x());
    int c_y = std::round(boundary_intercept.y());
    const int required_nodes = 10;
    int selected_nodes = 0;
    double selection_size = 0.0; //
    while (selected_nodes < required_nodes)
    {
        auto selection = std::pow(3.0 + 2.0 * selection_size, 2.0);
        int delta_range = std::static_cast<int> selection_size + 1.0;
        int span = delta_range * 2 + 1;


    }

}
