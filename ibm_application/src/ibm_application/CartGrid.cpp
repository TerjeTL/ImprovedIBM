
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

// method for calculating the pseudo-Inverse as recommended by Eigen developers
template<typename _Matrix_Type_>
_Matrix_Type_ PseudoInverseSVD(const _Matrix_Type_& a, double epsilon = std::numeric_limits<double>::epsilon())
{
    //Eigen::JacobiSVD< _Matrix_Type_ > svd(a, Eigen::ComputeFullU | Eigen::ComputeFullV);
    // For a non-square matrix
    Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = epsilon * std::max(a.cols(), a.rows()) * svd.singularValues().array().abs()(0);
    return svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}

double WeightingFunc(double x_n, double y_n, double a_d)
{
    double d_n = std::sqrt(std::pow(x_n, 2.0) + std::pow(y_n, 2.0));

    //std::cout << "relative dist: " << d_n << "\n\n";

    return std::exp(-std::pow(d_n, 2.0) / a_d);
}

void CartGrid::WLSQUpdateGeometry()
{
	for (auto& [ij, wlsq] : m_wlsq_data)
	{
        const int i = ij.first;
        const int j = ij.second;

        // get the selection of nodes

		// start by getting center location of node selection square (god gave us no solution for the circle problem)
        auto image_point = GetGridCoordinate(GetImagePoint(i, j));
        Eigen::Vector2d ghost_point = { static_cast<double>(i), static_cast<double>(j) };

        auto boundary_intercept = ghost_point + (image_point - ghost_point) / 2.0;
        auto ghost_pt_rel = ghost_point - boundary_intercept;

        // Find the dimensions needed to grab the required number of points
        int c_x = ghost_point.x();
        int c_y = ghost_point.y();
        const int required_nodes = 25;
        int selection_size = 0;

        while (wlsq.m_active_nodes_num < required_nodes)
        {
            // Define the size of the search area
            auto size = std::round(std::pow(3.0 + 2.0 * selection_size, 2.0));
            int delta_range = selection_size + 1;
            int span = delta_range * 2 + 1;

            wlsq.m_x_corner = std::clamp(c_x - delta_range, 0, (int)grid_flags.cols() - span);
            wlsq.m_y_corner = std::clamp(c_y - delta_range, 0, (int)grid_flags.rows() - span);

            // Get the subselection of nodes and check number of fluid points
            wlsq.m_subgrid = grid_flags.block(wlsq.m_y_corner, wlsq.m_x_corner, span, span);
            //std::cout << selection << "\n\n";
            wlsq.m_active_nodes_num = (wlsq.m_subgrid.array() == 0).count();

            selection_size++;
        }

        // construct vandermonde matrix, starting with GP as first row
        wlsq.m_vandermonde = Eigen::MatrixXd::Zero(wlsq.m_active_nodes_num + 1, 10);

        //vandermonde.row(0) = Eigen::Vector4d{ 1.0, 0.0, 0.0 };
        wlsq.m_vandermonde.row(0) = Eigen::Vector<double, 10>{
            1.0, ghost_pt_rel.x(), ghost_pt_rel.y(), ghost_pt_rel.x() * ghost_pt_rel.y(),
            std::pow(ghost_pt_rel.x(), 2.0), std::pow(ghost_pt_rel.y(), 2.0), std::pow(ghost_pt_rel.x(), 2.0) * ghost_pt_rel.y(), std::pow(ghost_pt_rel.y(), 2.0) * ghost_pt_rel.x(),
            std::pow(ghost_pt_rel.x(), 3.0), std::pow(ghost_pt_rel.y(), 3.0)
        };

        // need to get x's and y's relative to the boundary intercept for each selected fluid node
        int idx = 1;
        for (size_t i = 0; i < wlsq.m_subgrid.cols(); i++)
        {
            for (size_t j = 0; j < wlsq.m_subgrid.rows(); j++)
            {
                if (wlsq.m_subgrid(j, i) == 0)
                {
                    Eigen::Vector2d selected_node_rel_loc = Eigen::Vector2d{ static_cast<double>(wlsq.m_x_corner + i), static_cast<double>(wlsq.m_y_corner + j) } - boundary_intercept;
                    wlsq.m_vandermonde.row(idx) = Eigen::Vector<double, 10>{
                        1.0, selected_node_rel_loc.x(), selected_node_rel_loc.y(), selected_node_rel_loc.x() * selected_node_rel_loc.y(),
                        std::pow(selected_node_rel_loc.x(), 2.0), std::pow(selected_node_rel_loc.y(), 2.0), std::pow(selected_node_rel_loc.x(), 2.0) * selected_node_rel_loc.y(), std::pow(selected_node_rel_loc.y(), 2.0) * selected_node_rel_loc.x(),
                        std::pow(selected_node_rel_loc.x(), 3.0), std::pow(selected_node_rel_loc.y(), 3.0)
                    };

                    idx++;
                }
            }
        }

        // Construct weighting matrix
        Eigen::VectorXd weighting_coeffs(wlsq.m_active_nodes_num + 1);

        wlsq.m_a_d = std::pow(ghost_pt_rel.x(), 2.0) + std::pow(ghost_pt_rel.y(), 2.0);
        for (size_t i = 0; i < wlsq.m_subgrid.cols(); i++)
        {
            for (size_t j = 0; j < wlsq.m_subgrid.rows(); j++)
            {
                if (wlsq.m_subgrid(j, i) == 0)
                {
                    Eigen::Vector2d selected_node_rel_loc = Eigen::Vector2d{ static_cast<double>(wlsq.m_x_corner + i), static_cast<double>(wlsq.m_y_corner + j) } - boundary_intercept;

                    wlsq.m_a_d += std::pow(selected_node_rel_loc.x(), 2.0) + std::pow(selected_node_rel_loc.y(), 2.0);
                }
            }
        }
        wlsq.m_a_d = wlsq.m_a_d * m_weight_scaling; /// static_cast<double>(wlsq.m_active_nodes_num * wlsq.m_weight_scaling);

        wlsq.dist.clear();
        wlsq.dist.reserve(wlsq.m_active_nodes_num);
        wlsq.weight.clear();
        wlsq.weight.reserve(wlsq.m_active_nodes_num);

        weighting_coeffs[0] = WeightingFunc(ghost_pt_rel.x(), ghost_pt_rel.y(), wlsq.m_a_d);
        //weighting_coeffs[1] = WeightingFunc(ghost_pt_rel.x(), ghost_pt_rel.y(), a_d);

        int n = 1;
        for (size_t i = 0; i < wlsq.m_subgrid.cols(); i++)
        {
            for (size_t j = 0; j < wlsq.m_subgrid.rows(); j++)
            {
                if (wlsq.m_subgrid(j, i) == 0)
                {
                    Eigen::Vector2d selected_node_rel_loc = Eigen::Vector2d{ static_cast<double>(wlsq.m_x_corner + i), static_cast<double>(wlsq.m_y_corner + j) } - boundary_intercept;
                    weighting_coeffs[n] = WeightingFunc(selected_node_rel_loc.x(), selected_node_rel_loc.y(), wlsq.m_a_d);

                    // debug
                    wlsq.dist.push_back(std::sqrt(std::pow(selected_node_rel_loc.x(), 2.0) + std::pow(selected_node_rel_loc.y(), 2.0)));
                    wlsq.weight.push_back(weighting_coeffs[n]);

                    n++;
                }
            }
        }

        wlsq.m_weight = weighting_coeffs.asDiagonal();

        Eigen::MatrixXd w_v_product = wlsq.m_weight * wlsq.m_vandermonde;
        //wlsq.m_M = w_v_product.completeOrthogonalDecomposition().pseudoInverse() * wlsq.m_weight;

        wlsq.m_M = PseudoInverseSVD(w_v_product) * wlsq.m_weight;
	}
}

void CartGrid::WeightedLeastSquaresMethod()
{
    m_wlsq_phi_matrix = phi_matrix;

    for (auto& [ij, wlsq] : m_wlsq_data)
    {
        int i = ij.first;
        int j = ij.second;

        // construct phi vector
        Eigen::VectorXd phi_vec(wlsq.m_active_nodes_num + 1);

        size_t parent_sdf = ghost_point_parent_sdf(j, i);
        phi_vec[0] = m_wlsq_phi_matrix(j, i);

        int n = 1;
        for (size_t i = 0; i < wlsq.m_subgrid.cols(); i++)
        {
            for (size_t j = 0; j < wlsq.m_subgrid.rows(); j++)
            {
                if (wlsq.m_subgrid(j, i) == 0)
                {
                    phi_vec[n] = m_wlsq_phi_matrix(wlsq.m_y_corner + j, wlsq.m_x_corner + i);
                    n++;
                }
            }
        }

        // solve system Ax = b
        Eigen::VectorXd c_vec = wlsq.m_M * phi_vec;

        double test = 0.0;
        if (GetBoundaryCondition(i, j) == BoundaryCondition::Dirichlet)
        {
            for (size_t n = 1; n < wlsq.m_M.cols(); n++)
            {
                Eigen::VectorXd pos = wlsq.m_vandermonde.row(n);
                double phi_int = c_vec.dot(pos);
                test += wlsq.m_M(0, n) * phi_int;
            }

            wlsq.m_gp_val = (immersed_boundaries.at(parent_sdf)->GetBoundaryPhi() - test) / wlsq.m_M(0, 0);
        }
        else if (GetBoundaryCondition(i, j) == BoundaryCondition::Neumann)
        {
            auto image_point = GetGridCoordinate(GetImagePoint(i, j));
            Eigen::Vector2d ghost_point = { static_cast<double>(i), static_cast<double>(j) };

            Eigen::Vector2d normal = (image_point - ghost_point).normalized();
            auto n_x = normal.x();
            auto n_y = normal.y();

            //std::cout << "(i, j): " << i << ", " << j << "\n\n";
            //std::cout << "normal: " << n_x << ", " << n_y << "\n\n";

            for (size_t n = 1; n < wlsq.m_M.cols(); n++)
            {
                Eigen::VectorXd pos = wlsq.m_vandermonde.row(n);
                double phi_int = c_vec.dot(pos);
                test += (n_x * wlsq.m_M(1, n) + n_y * wlsq.m_M(2, n)) * phi_int;
            }

            //immersed_boundaries.at(parent_sdf)->GetNormal()

            wlsq.m_gp_val = (immersed_boundaries.at(parent_sdf)->GetBoundaryPhi()/(double)(GetMeshSize().first - 1) - test) / (n_x * wlsq.m_M(1, 0) + n_y * wlsq.m_M(2, 0));
            /*std::cout << "PHI VECTOR\n" << phi_vec << "\n\n";
            std::cout << "M\n" << wlsq.m_M << "\n\n";
            std::cout << "COEFFICIENT VEC\n" << c_vec << "\n\n";*/
        }

        
        // A = V^T W V
        //auto A = vandermonde.transpose() * (weighting_matrix) * vandermonde;
        //std::cout << "A\n" << A << "\n\n";

        // b = V^T W Phi
        //auto b = vandermonde.transpose() * (weighting_matrix) * phi_vec;
        //std::cout << "b\n" << b << "\n\n";

        //Eigen::Vector4d C = A.colPivHouseholderQr().solve(b);
        // SVD: Eigen::Vector4d C = (weighting_matrix * vandermonde).template bdcSvd<Eigen::ComputeThinU | Eigen::ComputeThinV>().solve(weighting_matrix * phi_vec);
        //std::cout << "C\n" << C << "\n\n";

        //Eigen::Vector4d gp_vec{ 1.0, ghost_pt_rel.x(), ghost_pt_rel.y(), ghost_pt_rel.x() * ghost_pt_rel.y() };

#ifdef DEBUG_TERMINAL_WLSQ
        std::cout << "IP\n" << image_point << "\n\n";
        std::cout << "GP\n" << ghost_point << "\n\n";
        std::cout << "BI\n" << boundary_intercept << "\n\n";
        std::cout << "FLAG MATRIX\n" << grid_flags << "\n\n";
        std::cout << "SELECTION\n" << selection << "\n\n";
        std::cout << "VANDERMONDE\n" << vandermonde << "\n\n";
        std::cout << "WEIGHTING MATRIX\n" << weighting_matrix << "\n\n";
        std::cout << "PHI VECTOR\n" << phi_vec << "\n\n";
        std::cout << "M\n" << M << "\n\n";
        std::cout << "Mt\n" << M.transpose() << "\n\n";
        std::cout << "COEFFICIENT VEC\n" << c_vec << "\n\n";
        //std::cout << "GP VECTOR\n" << gp_vec << "\n\n";
        std::cout << "GP RECONSTRUCTION\n" << gp_val << "\n\n";
#endif
    }
}
