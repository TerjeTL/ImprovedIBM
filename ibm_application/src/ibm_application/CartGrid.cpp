
// This file defines the cartesian grid

#include "ibm_application/CartGrid.h"

#include<Eigen/Core>
#include<Eigen/SVD>

#include <iostream>
#include <assert.h>

//#define DEBUG_TERMINAL_WLSQ

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

    WLSQInit();
}


// Apply linear distribution to the initialization
void CartGrid::InitializeFieldUnsteady(std::vector<double>& phi_analytical)
{
    assert(immersed_boundaries.size() == 1 && "Field initialization does not match case setup");

    auto r_outer = immersed_boundaries.begin()->second->GetSize();
    const Eigen::MatrixXi& grid_flags = GetGridFlags();
    Eigen::MatrixXd analytical = Eigen::MatrixXd::Zero(grid_flags.rows(), grid_flags.cols());

    for (int j = 0; j < grid_flags.rows(); ++j)
    {
        for (int i = 0; i < grid_flags.cols(); ++i)
        {
            auto coordinate = GetWorldCoordinate({ i, j });

            if (grid_flags(j, i) == 0)
            {
                auto r = std::sqrt(std::pow(coordinate.x() - 0.5, 2) + std::pow(coordinate.y() - 0.5, 2));

                auto dr = r_outer / (phi_analytical.size() - 1);
                int idx_min = static_cast<int>(r / dr);

                double phi_interp = std::lerp(phi_analytical[idx_min], phi_analytical[idx_min + 1], r - idx_min * dr);

                analytical(j, i) = phi_interp * 2.0;
            }
        }
    }

    phi_matrix = analytical;
}

// Apply linear distribution to the initialization
void CartGrid::InitializeField()
{
    assert(immersed_boundaries.size() == 2 && "Can only do custom field initialization with exactly two immersed boundaries");


    // Steady state case
    for (size_t j = 0UL; j < grid_flags.rows(); ++j) {
        for (size_t i = 0UL; i < grid_flags.cols(); ++i) {

            // Skip inactive cells 
            if ( grid_flags(j, i) == 1 )
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


// Implementing WLSQ method
// NB! Make sure to call init before using
void CartGrid::WLSQInit()
{
    m_wlsq_phi_matrix = phi_matrix;

    if (m_wlsq_data.empty())
    {
	    for (size_t i = 0; i < grid_flags.cols(); ++i)
	    {
		    for (size_t j = 0; j < grid_flags.rows(); ++j)
		    {
			    if (grid_flags(j,i) == 2)
			    {
                    m_wlsq_data[{i, j}] = WLSQdata{};
			    }		    
		    }
	    }
    }

    WLSQUpdateGeometry();
}

void CartGrid::ConstructNumericalStencil(WLSQdata& data, size_t required_nodes)
{
    auto image_point = GetGridCoordinate(GetImagePoint(data.ghost_point.x(), data.ghost_point.y()));
    Eigen::Vector2d boundary_intercept_grid_coord = data.ghost_point.cast<double>() + (image_point - data.ghost_point.cast<double>()) / 2.0;

    // Find the dimensions needed to grab the required number of points
    int c_x = (int)std::round(boundary_intercept_grid_coord.x());
    int c_y = (int)std::round(boundary_intercept_grid_coord.y());

    int selection_size = 0;

    int x_corner = 0;
    int y_corner = 0;

    while (data.m_active_nodes_num < required_nodes)
    {
        // Define the size of the search area
        auto size = std::round(std::pow(3.0 + 2.0 * selection_size, 2.0));
        int delta_range = selection_size + 1;
        int span = delta_range * 2 + 1;

        x_corner = std::clamp(c_x - delta_range, 0, (int)grid_flags.cols() - span);
        y_corner = std::clamp(c_y - delta_range, 0, (int)grid_flags.rows() - span);

        // Get the subselection of nodes and check number of fluid points
        data.m_subgrid = grid_flags.block(y_corner, x_corner, span, span);
        //std::cout << selection << "\n\n";
        data.m_active_nodes_num = (data.m_subgrid.array() == 0).count();

        selection_size++;
    }

    auto body_intercpt_wrld = GetWorldCoordinate(boundary_intercept_grid_coord);
    auto wrld_loc = GetWorldCoordinate( {data.ghost_point.x(), data.ghost_point.y() });

    data.m_pos.clear();
    data.m_pos.reserve(data.m_active_nodes_num + 1);
    data.m_pos.push_back(wrld_loc - body_intercpt_wrld);

    data.m_numerical_stencil.clear();
    data.m_numerical_stencil.reserve(data.m_active_nodes_num + 1);
    data.m_numerical_stencil.push_back({ data.ghost_point.x(), data.ghost_point.y() });

    for (size_t i = 0; i < data.m_subgrid.cols(); i++) {
        for (size_t j = 0; j < data.m_subgrid.rows(); j++) {
            if (data.m_subgrid(j, i) == 0) {
                auto wrld_pos = GetWorldCoordinate(Eigen::Vector2d{ static_cast<double>(x_corner + i), static_cast<double>(y_corner + j) });
                auto wrld_pos_rel = wrld_pos - body_intercpt_wrld;

                data.m_pos.push_back(wrld_pos_rel);
                data.m_numerical_stencil.push_back({ x_corner + i, y_corner + j });
            }
        }
    }
}

// method for calculating the pseudo-Inverse as recommended by Eigen developers
template<typename _Matrix_Type_>
std::pair<_Matrix_Type_, double> PseudoInverseSVD(const _Matrix_Type_& a, double epsilon = std::numeric_limits<double>::epsilon())
{
    //Eigen::JacobiSVD< _Matrix_Type_ > svd(a, Eigen::ComputeFullU | Eigen::ComputeFullV);
    // For a non-square matrix
    Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = epsilon * std::max(a.cols(), a.rows()) * svd.singularValues().array().abs()(0);

    double condition_number = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);

    return { svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint(), condition_number };
}

double WeightingFunc(double x_n, double y_n, double a_d)
{
    double d_n = std::sqrt(std::pow(x_n, 2.0) + std::pow(y_n, 2.0));

    //std::cout << "relative dist: " << d_n << "\n\n";

    return std::exp(-std::pow(d_n, 2.0) / a_d);
}

Eigen::MatrixXd ConstructVandermonde(int polynomial_order, CartGrid::WLSQdata& state, const std::vector<Eigen::Vector2d>& relative_positions)
{
    Eigen::MatrixXd vandermonde;

	switch (polynomial_order)
	{
	case 1:
	{
		vandermonde = Eigen::MatrixXd::Zero(state.m_active_nodes_num + 1, 3);

        for (int n = 0; n < relative_positions.size(); ++n)
        {
            const auto& rel_pos = relative_positions[n];

            vandermonde.row(n) = Eigen::Vector<double, 3>{
                        1.0, rel_pos.x(), rel_pos.y()
            };
        }
    break;
	}
	case 2:
    {
        vandermonde = Eigen::MatrixXd::Zero(state.m_active_nodes_num + 1, 6);

        for (int n = 0; n < relative_positions.size(); ++n)
        {
            const auto& rel_pos = relative_positions[n];

            vandermonde.row(n) = Eigen::Vector<double, 6>{
                        1.0, rel_pos.x(), rel_pos.y(), rel_pos.x() * rel_pos.y(),
                        std::pow(rel_pos.x(), 2.0), std::pow(rel_pos.y(), 2.0)
            };
        }
        break;
    }
    case 3:
    {
        vandermonde = Eigen::MatrixXd::Zero(state.m_active_nodes_num + 1, 10);

        for (int n = 0; n < relative_positions.size(); ++n)
        {
            const auto& rel_pos = relative_positions[n];

            vandermonde.row(n) = Eigen::Vector<double, 10>{
                    1.0, rel_pos.x(), rel_pos.y(), rel_pos.x() * rel_pos.y(),
                    std::pow(rel_pos.x(), 2.0), std::pow(rel_pos.y(), 2.0),
                    std::pow(rel_pos.x(), 2.0) * rel_pos.y(), std::pow(rel_pos.y(), 2.0) * rel_pos.x(),
                    std::pow(rel_pos.x(), 3.0), std::pow(rel_pos.y(), 3.0)
            };
        }
        break;
    }
	default:
        assert("Invalid Polynomial Order");
	}

    return vandermonde;
}

Eigen::MatrixXd ConstructWeightMatrix(double weight_scaling, CartGrid::WLSQdata& state, const std::vector<Eigen::Vector2d>& relative_positions)
{
    // Construct weighting matrix
    Eigen::VectorXd weighting_coeffs(state.m_active_nodes_num + 1);

    double a_d = 0.0;
    for (auto& rel_pos : relative_positions)
    {
        a_d += std::pow(rel_pos.x(), 2.0) + std::pow(rel_pos.y(), 2.0);
    }
    a_d *= weight_scaling;

    state.dist.clear();
    state.dist.reserve(state.m_active_nodes_num);
    state.weight.clear();
    state.weight.reserve(state.m_active_nodes_num);

    for (int n = 0; n < relative_positions.size(); ++n)
    {
        weighting_coeffs[n] = WeightingFunc(relative_positions[n].x(), relative_positions[n].y(), a_d);

        if (n > 0)
        {
            state.dist.push_back(std::sqrt(std::pow(relative_positions[n].x(), 2.0) + std::pow(relative_positions[n].y(), 2.0)));
            state.weight.push_back(weighting_coeffs[n]);
        }
    }

    return weighting_coeffs.asDiagonal();
}

void CartGrid::WLSQUpdateGeometry()
{
    double cond_max = 0.0;

	for (auto& [ij, wlsq] : m_wlsq_data)
	{
        const int i = ij.first;
        const int j = ij.second;

        wlsq.ghost_point = { i, j };

        size_t parent_sdf = ghost_point_parent_sdf(j, i);
        Eigen::Vector2d world_loc = GetWorldCoordinate(Eigen::Vector2d{ i, j });
        Eigen::Vector2d unit_normal = immersed_boundaries.at(parent_sdf)->GetNormal(world_loc.x(), world_loc.y());

        // get the selection of nodes

		// start by getting center location of node selection square (god gave us no solution for the circle problem)
        const int required_nodes = 25;


        ConstructNumericalStencil(wlsq, required_nodes);

        // construct vandermonde matrix, starting with GP as first row
        wlsq.m_vandermonde = ConstructVandermonde(3, wlsq, wlsq.m_pos);

        wlsq.m_weight = ConstructWeightMatrix(m_weight_scaling, wlsq, wlsq.m_pos);

        Eigen::MatrixXd w_v_product = wlsq.m_weight * wlsq.m_vandermonde;

        auto pseudo_inv = PseudoInverseSVD(w_v_product);
        cond_max = std::max(cond_max, pseudo_inv.second);

        wlsq.m_M = pseudo_inv.first * wlsq.m_weight;

        // Initialize related vars
        wlsq.m_phi_vec = Eigen::VectorXd::Zero(wlsq.m_active_nodes_num + 1);
        wlsq.m_bc_type = GetBoundaryCondition(i, j);
        wlsq.m_unit_normal = unit_normal;
        wlsq.m_bc_value = GetBoundaryPhi(i, j);

        // optimization measures (pre-slicing data)
        wlsq.m_num_stencil_reduced = std::vector<std::pair<int, int>>(wlsq.m_numerical_stencil.begin() + 1, wlsq.m_numerical_stencil.end());
        wlsq.m_phi_vec_reduced = wlsq.m_phi_vec(Eigen::seq(1, Eigen::placeholders::last));
        wlsq.m_M_0 = wlsq.m_M(0, Eigen::seq(1, Eigen::placeholders::last));
        wlsq.m_M_1 = wlsq.m_M(1, Eigen::seq(1, Eigen::placeholders::last));
        wlsq.m_M_2 = wlsq.m_M(2, Eigen::seq(1, Eigen::placeholders::last));

        if (wlsq.m_bc_type == BoundaryCondition::Dirichlet)
        {
            wlsq.m_M_den = 1.0 / wlsq.m_M(0, 0);
            wlsq.m_bc_term = wlsq.m_bc_value * wlsq.m_M_den;
        }
        else if (wlsq.m_bc_type == BoundaryCondition::Neumann)
        {
            //wlsq.m_M_den = 1.0 / wlsq.m_M(0, 0);
            //wlsq.m_bc_term = wlsq.m_bc_value * wlsq.m_M_den;
            //(n_x * wlsq.m_M(1, 0) + n_y * wlsq.m_M(2, 0));
        }
	}

    printf("Condition Number: %g\n\n", cond_max);
}

//-------------------
// WLSQ Main Update
//-------------------
// This function updates all ghost point values using the WLSQ method. Prerequisites for
// this method is that the wlsq data of each ghost point has been previously initialized
// with WLSQInit, and that WLSQUpdateGeometry has been called if there have been any changes
// to the configuration such as the weight function or boundary geometries.
void CartGrid::WeightedLeastSquaresMethod()
{
    
    for (auto& [ij, wlsq] : m_wlsq_data)
    {
        // Update vector of phi values corresponding to the node selection
        int n = 0;
		for (auto [i, j] : wlsq.m_num_stencil_reduced) // pre-sliced numerical stencil
		{
            wlsq.m_phi_vec_reduced[n] = phi_matrix(j, i);
            n++;
		}

        // Solve based on boundary condition
        // - See pp. ...
        double linear_comb_sum = 0.0;
        if (wlsq.m_bc_type == BoundaryCondition::Dirichlet)
        {
            //for (size_t n = 1; n < wlsq.m_M.cols(); n++)
            //{
            //    linear_comb_sum += wlsq.m_M(0, n) * wlsq.phi_vec(n);
            //}
            linear_comb_sum = wlsq.m_M_0.dot(wlsq.m_phi_vec_reduced);

            //wlsq.m_gp_val = (wlsq.m_bc_value - linear_comb_sum) / wlsq.m_M(0, 0);
            wlsq.m_gp_val = wlsq.m_bc_term - linear_comb_sum * wlsq.m_M_den;
        }
        else if (wlsq.m_bc_type == BoundaryCondition::Neumann)
        {
            auto n_x = wlsq.m_unit_normal.x();
            auto n_y = wlsq.m_unit_normal.y();

            for (size_t n = 1; n < wlsq.m_M.cols(); n++)
            {
                linear_comb_sum += (n_x * wlsq.m_M(1, n) + n_y * wlsq.m_M(2, n)) * wlsq.m_phi_vec(n);
            }
            //double test_linear_comb_sum = (n_x * wlsq.m_M.row(1)(Eigen::seq(1, Eigen::placeholders::last)) + n_y * wlsq.m_M.row(2)(Eigen::seq(1, Eigen::placeholders::last))) * wlsq.phi_vec(Eigen::seq(1, Eigen::placeholders::last));

            wlsq.m_gp_val = (wlsq.m_bc_value - linear_comb_sum) / (n_x * wlsq.m_M(1, 0) + n_y * wlsq.m_M(2, 0));
        }
    }
}
