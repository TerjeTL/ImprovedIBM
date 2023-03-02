
// This file defines the cartesian grid

#include "ibm_application/CartGrid.h"

#include <Eigen/Dense>

#include <iostream>
#include <assert.h>

//#define DEBUG_TERMINAL

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
            for (size_t j = 0UL; j < grid_flags.rows(); ++j) {
                for (size_t i = 0UL; i < grid_flags.cols(); ++i) {

                    auto x = static_cast<double>(i) * h;
                    auto y = static_cast<double>(j) * h;

                    // Check if i,j should be a ghost, active or inactive grid point and assign correct flag
                    if (sdf->SignedDistanceFunction(x, y) < 0.0) {
                        if (sdf->SignedDistanceFunction(x + h, y) >= 0.0
                            || sdf->SignedDistanceFunction(x - h, y) >= 0.0
                            || sdf->SignedDistanceFunction(x, y + h) >= 0.0
                            || sdf->SignedDistanceFunction(x, y - h) >= 0.0) {
                            
                            // Set location descriptors
                            grid_flags(j, i) = 2;
                            ghost_point_parent_sdf(j, i) = hash;

                            // Calculate normal and assign respective image-point (they are not cleared!)
                            auto normal = sdf->GetNormal(x, y);
                            Eigen::Vector2d node{ x, y };

                            Eigen::Vector2d image_point = node + normal * std::abs(sdf->SignedDistanceFunction(x, y)) * m_ip_stencil_length_factor;
                            image_points[{ i, j }] = image_point;

                            // Regenerate coefficients for bilinear interpolation
                            UpdateInterpolationCoeffs(i, j);

                            phi_image_point_matrix(j, i) = BilinearInterpolation(i, j);
                            boundary_phi(j, i) = sdf->GetBoundaryPhi();
                        }
                        else
                        {
                            // Set cell as inactive (solid)
                            grid_flags(j, i) = 1;
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

    for (size_t j = 0UL; j < grid_flags.rows(); ++j) {
        for (size_t i = 0UL; i < grid_flags.cols(); ++i) {

            // Skip inactive cells 
            if (grid_flags(j, i) != 0)
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
            phi_matrix(j,i) = dist_1 / std::abs(dist_1 + dist_2) * (phi_2 - phi_1) + phi_1;
        }
    }
}

// Update beta coefficients
void CartGrid::UpdateInterpolationCoeffs(size_t i, size_t j)
{
    // Interpolation in grid-space
    Eigen::Vector2d world_coordinate = image_points.at({ i, j });
    Eigen::Vector2d grid_coordinate = GetGridCoordinate(image_points.at({ i, j }));
    size_t parent_sdf = ghost_point_parent_sdf(j, i);
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

        if (grid_flags(node_indices[p](1), node_indices[p](0)) == 2)
        {
            auto image_pt_loc = GetGridCoordinate(GetImagePoint(node_indices[p](0), node_indices[p](1)));

            auto world_normal = immersed_boundaries.at(parent_sdf)->GetNormal(world_coordinate.x(), world_coordinate.y());
            world_normal.normalize();
            auto grid_normal = GetGridCoordinate(world_normal);

            node_locations[p] = node_indices[p].cast<double>() + (image_pt_loc - node_indices[p].cast<double>()) / 2.0;

            if (boundary == BoundaryCondition::Dirichlet)
            {
                node_phi(p) = boundary_phi(node_indices[p](1), node_indices[p](0));
                vandermonde_row = { 1, node_locations[p](0), node_locations[p](1), node_locations[p](0) * node_locations[p](1) };
            }
            else // Neumann
            {
                node_phi(p) = boundary_phi(node_indices[p](1), node_indices[p](0));
                vandermonde_row = { 0, grid_normal.x(), grid_normal.y(), node_locations[p](0) * grid_normal.y() + node_locations[p](1) * grid_normal.x() };
            }
        }
        else
        {
            node_locations[p] = node_indices[p].cast<double>();
            node_phi(p) = phi_matrix(node_indices[p](1), node_indices[p](0));
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
    size_t parent_sdf = ghost_point_parent_sdf(j, i);
    BoundaryCondition boundary = immersed_boundaries.at(parent_sdf)->GetBoundaryCondition();

    std::array<Eigen::Vector2i, 4> node_indices;
    Eigen::Vector4d node_phi;

    node_indices[0] = grid_coordinate.cast<int>();
    node_indices[1] = node_indices[0] + Eigen::Vector2i{ 0,1 };
    node_indices[2] = node_indices[0] + Eigen::Vector2i{ 1,1 };
    node_indices[3] = node_indices[0] + Eigen::Vector2i{ 1,0 };

    for (size_t p = 0; p < node_indices.size(); p++)
    {
        if (grid_flags(node_indices[p](1), node_indices[p](0)) == 2)
        {
            auto image_pt_loc = GetGridCoordinate(GetImagePoint(node_indices[p](0), node_indices[p](1)));
            if (boundary == BoundaryCondition::Dirichlet)
            {
                node_phi(p) = boundary_phi(node_indices[p](1), node_indices[p](0));
            }
            else
            {
                node_phi(p) = boundary_phi(node_indices[p](1), node_indices[p](0));
            }
        }
        else
        {
            node_phi(p) = phi_matrix(node_indices[p](1), node_indices[p](0));
        }
    }

    Eigen::Vector4d coefficients_new = beta_coeffs.at({ i,j });
    double phi = coefficients_new.dot(node_phi);

    return phi;
}


double WeightingFunc(double x_n, double y_n, double a_d)
{
    double d_n = std::sqrt(std::pow(x_n, 2.0) + std::pow(y_n, 2.0));

    //std::cout << "relative dist: " << d_n << "\n\n";

    return std::exp(-std::pow(d_n, 2.0)/a_d);
}

// Implementing WLSQ method
double CartGrid::WeightedLeastSquaresMethod(size_t i, size_t j)
{
    // get the selection of nodes
    
    // start by getting center location of node selection square (god gave us no solution for the circle problem)
    auto image_point = GetGridCoordinate(GetImagePoint(i, j));
    Eigen::Vector2d ghost_point = { static_cast<double>(i), static_cast<double>(j) };

    auto boundary_intercept = ghost_point + (image_point - ghost_point) / 2.0;
    auto ghost_pt_rel = ghost_point - boundary_intercept;

    // Find the dimensions needed to grab the required number of points
    int c_x = ghost_point.x();
    int c_y = ghost_point.y();
    const int required_nodes = 15;
    Eigen::MatrixXi selection;
    int selected_nodes = 0;
    int selection_size = 0;
    int x_corner = 0;
    int y_corner = 0;

    while (selected_nodes < required_nodes)
    {
        // Define the size of the search area
        auto size = std::round(std::pow(3.0 + 2.0 * selection_size, 2.0));
        int delta_range = selection_size + 1;
        int span = delta_range * 2 + 1;

        x_corner = std::clamp(c_x - delta_range, 0, (int)grid_flags.cols() - span - 1 );
        y_corner = std::clamp(c_y - delta_range, 0, (int)grid_flags.rows() - span - 1 );

        // Get the subselection of nodes and check number of fluid points
        selection = grid_flags.block(y_corner, x_corner, span, span);
        //std::cout << selection << "\n\n";
        selected_nodes = (selection.array() == 0).count();

        selection_size++;
    }

    // construct vandermonde matrix, starting with BI as first row, GP as second row
    Eigen::MatrixXd vandermonde = Eigen::MatrixXd::Zero(selected_nodes + 1, 4);
    
    vandermonde.row(0) = Eigen::Vector4d{ 1.0, 0.0, 0.0, 0.0 };
    //vandermonde.row(1) = Eigen::Vector4d{ 1.0, GetWorldCoordinate(ghost_pt_rel).x(), GetWorldCoordinate(ghost_pt_rel).y(), GetWorldCoordinate(ghost_pt_rel).x() * GetWorldCoordinate(ghost_pt_rel).y() };

    // need to get x's and y's relative to the boundary intercept for each selected fluid node
    int idx = 1;
    for (size_t i = 0; i < selection.cols(); i++)
    {
        for (size_t j = 0; j < selection.rows(); j++)
        {
            if (selection(j, i) == 0)
            {
                Eigen::Vector2d selected_node_rel_loc = Eigen::Vector2d{ static_cast<double>(x_corner + i), static_cast<double>(y_corner + j) } - boundary_intercept;
                vandermonde.row(idx) = Eigen::Vector4d{ 1.0, selected_node_rel_loc.x(), selected_node_rel_loc.y(), selected_node_rel_loc.x() * selected_node_rel_loc.y() };

                idx++;
            }
        }
    }

    // Construct weighting matrix
    Eigen::VectorXd weighting_coeffs(selected_nodes + 1);

    double a_d = 0.0; // std::pow(GetWorldCoordinate(ghost_pt_rel).x(), 2.0) + std::pow(GetWorldCoordinate(ghost_pt_rel).y(), 2.0);
    for (size_t i = 0; i < selection.cols(); i++)
    {
        for (size_t j = 0; j < selection.rows(); j++)
        {
            if (selection(j, i) == 0)
            {
                Eigen::Vector2d selected_node_rel_loc = Eigen::Vector2d{ static_cast<double>(x_corner + i), static_cast<double>(y_corner + j) } - boundary_intercept;

                a_d += std::pow(selected_node_rel_loc.x(), 2.0) + std::pow(selected_node_rel_loc.y(), 2.0);
            }
        }
    }
    a_d = a_d / static_cast<double>(selected_nodes);

    weighting_coeffs[0] = 1.0;
    //weighting_coeffs[1] = WeightingFunc(ghost_pt_rel.x(), ghost_pt_rel.y(), a_d);

    int n = 1;
    for (size_t i = 0; i < selection.cols(); i++)
    {
        for (size_t j = 0; j < selection.rows(); j++)
        {
            if (selection(j, i) == 0)
            {
                Eigen::Vector2d selected_node_rel_loc = Eigen::Vector2d{ static_cast<double>(x_corner + i), static_cast<double>(y_corner + j) } - boundary_intercept;
                weighting_coeffs[n] = WeightingFunc(selected_node_rel_loc.x(), selected_node_rel_loc.y(), a_d);

                n++;
            }
        }
    }
    
    Eigen::MatrixXd weighting_matrix = weighting_coeffs.asDiagonal();

    // construct phi vector
    Eigen::VectorXd phi_vec(selected_nodes + 1);

    size_t parent_sdf = ghost_point_parent_sdf(j, i);
    BoundaryCondition boundary = immersed_boundaries.at(parent_sdf)->GetBoundaryCondition();
    phi_vec[0] = immersed_boundaries.at(parent_sdf)->GetBoundaryPhi();
    //phi_vec[1] = phi_matrix(j, i);

    n = 1;
    for (size_t i = 0; i < selection.cols(); i++)
    {
        for (size_t j = 0; j < selection.rows(); j++)
        {
            if (selection(j, i) == 0)
            {
                phi_vec[n] = phi_matrix(y_corner + j, x_corner + i);
                n++;
            }
        }
    }

    // solve system Ax = b

    // A = V^T W V
    auto A = vandermonde.transpose() * (weighting_matrix) * vandermonde;
    //std::cout << "A\n" << A << "\n\n";

    // b = V^T W Phi
    auto b = vandermonde.transpose() * (weighting_matrix) * phi_vec;
    //std::cout << "b\n" << b << "\n\n";

    Eigen::Vector4d C = A.colPivHouseholderQr().solve(b);
    // SVD: Eigen::Vector4d C = (weighting_matrix * vandermonde).template bdcSvd<Eigen::ComputeThinU | Eigen::ComputeThinV>().solve(weighting_matrix * phi_vec);
    //std::cout << "C\n" << C << "\n\n";

    Eigen::Vector4d gp_vec{ 1.0, ghost_pt_rel.x(), ghost_pt_rel.y(), ghost_pt_rel.x() * ghost_pt_rel.y() };

    auto res = C.dot(gp_vec);

#ifdef DEBUG_TERMINAL
    std::cout << "IP\n" << image_point << "\n\n";
    std::cout << "GP\n" << ghost_point << "\n\n";
    std::cout << "BI\n" << boundary_intercept << "\n\n";
    std::cout << "FLAG MATRIX\n" << grid_flags << "\n\n";
    std::cout << "SELECTION\n" << selection << "\n\n";
    std::cout << "VANDERMONDE\n" << vandermonde << "\n\n";
    std::cout << "WEIGHTING MATRIX\n" << weighting_matrix << "\n\n";
    std::cout << "PHI VECTOR\n" << phi_vec << "\n\n";
    std::cout << "GP VECTOR\n" << gp_vec << "\n\n";
    std::cout << "GP RECONSTRUCTION\n" << res << "\n\n";

    auto normal = ghost_pt_rel.normalized();
    std::cout << "NORMAL DIR\n" << normal << "\n\n";
    std::vector<double> probe{ -1, -0.5, 0.0, 0.5, 1.0, 2.0 };
    std::vector<std::string> data;

    std::cout << "GP NORM POS\n" << ghost_pt_rel.dot(normal) << "\n\n";

    for (auto loc : probe)
    {
        auto pos = normal * loc;
        Eigen::Vector4d gp_vec{ 1.0, pos.x(), pos.y(), pos.x() * pos.y() };
        auto val = C.dot(gp_vec);

        auto str = "(" + std::to_string(loc) + "," + std::to_string(val) + ")";
        data.push_back(str);

        std::cout << str << "\n";
    }
#endif

    return res;
}
