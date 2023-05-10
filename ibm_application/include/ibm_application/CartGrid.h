#pragma once

#include <memory>
#include <functional>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "ibm_application/GeometrySDF.h"

class CartGrid
{
public:
	struct WLSQdata
	{
		double m_gp_val = 0.0;
		double m_a_d = 0.0;
		int m_active_nodes_num = 0;
		int m_x_corner = 0;
		int m_y_corner = 0;
		BoundaryCondition bc;

		std::vector<std::pair<int, int>> m_numerical_stencil; // Warning: ordered

		Eigen::Vector2d m_unit_normal;

		Eigen::MatrixXi m_subgrid;

		Eigen::MatrixXd m_vandermonde;
		Eigen::MatrixXd m_weight;
		Eigen::VectorXd phi_vec;
		Eigen::MatrixXd m_M;

		// debugging
		std::vector<double> dist;
		std::vector<double> weight;
	};

	CartGrid(size_t nn, double phi_0 = 0.0);
	~CartGrid() {};

	void InitializeFieldUnsteady(std::vector<double>& phi_analytical);
	void InitializeField();
	void UpdateGrid();

	template <typename T> int sgn(T val) {
		return (T(0) < val) - (val < T(0));
	}

	void AddImmersedBoundary(std::string name, std::shared_ptr<GeometrySDF> geometry)
	{
		if (geometry) {

			immersed_boundaries[std::hash<std::string>{}(name)] = geometry;
		}
	}

	int GetCellFlag(size_t i, size_t j) const
	{
		return grid_flags(j, i);
	}
	
	double GetPhi(size_t i, size_t j) const
	{
		return phi_matrix(j, i);
	}

	const std::vector<std::pair<int, int>> const GetGhostPoints()
	{
		std::vector<std::pair<int, int>> keys;
		keys.reserve(image_points.size());
		for (const auto& [key, val] : image_points)
		{
			keys.push_back(key);
		}
		return keys;
	}

	Eigen::Vector2d GetImagePoint(size_t i, size_t j) const
	{
		return image_points.at({ i, j });
	}
	double GetImagePointPhi(size_t i, size_t j) const
	{
		return phi_image_point_matrix(j, i);
	}
	double GetBoundaryPhi(size_t i, size_t j) const
	{
		return boundary_phi(j, i);
	}

	std::pair<int, int> GetMeshSize() const
	{
		return std::pair<int, int>{grid_flags.cols(), grid_flags.rows()};
	}

	std::unordered_map<std::size_t, std::shared_ptr<GeometrySDF>> GetImmersedBoundaries() const
	{
		return immersed_boundaries;
	}

	Eigen::Vector2d GetGridCoordinate(Eigen::Vector2d world_coordinate)
	{
		return world_coordinate / h;
	}

	Eigen::Vector2d GetWorldCoordinate(Eigen::Vector2d grid_coordinate) const
	{
		return grid_coordinate * h;
	}

	Eigen::Vector3d GetGridCellSize()
	{
		return length_scales * h;
	}

	void UpdateInterpolationCoeffs(size_t i, size_t j);

	double BilinearInterpolation(size_t i, size_t j);

	void WLSQInit();
	void WLSQUpdateGeometry();
	void WeightedLeastSquaresMethod();

	const Eigen::MatrixXd& GetPhiMatrix() const
	{
		return phi_matrix;
	}

	BoundaryCondition GetBoundaryCondition(size_t i, size_t j) const
	{
		return immersed_boundaries.at(ghost_point_parent_sdf(j, i))->GetBoundaryCondition();
	}

	Eigen::MatrixXd& GetPhiMatrixRef()
	{
		return phi_matrix;
	}

	Eigen::MatrixXi& GetGridFlags()
	{
		return grid_flags;
	}

	Eigen::MatrixXd& GetPhiImagePointMatrixRef()
	{
		return phi_image_point_matrix;
	}

	const std::map<std::pair<int, int>, WLSQdata>& GetWLSQdata() const { return m_wlsq_data; }
	
	// Debugging-oriented functions
	std::array<Eigen::Vector2d, 4> GetBilinearInterpSelection(size_t i, size_t j)
	{
		return bilinear_interp_selection.at({ i, j });
	}

	double GetPhiWLSQ(int i, int j)
	{
		return m_wlsq_data.at({ i, j }).m_gp_val;
	}

	static constexpr double m_ip_stencil_length_factor = 2.0; // determines how far into the domain the stencil goes. Based on GP-to-boundary length.

	double m_weight_scaling = 1.0;

private:
	std::unordered_map<std::size_t, std::shared_ptr<GeometrySDF>> immersed_boundaries;

	double h = 0.0;

	Eigen::Vector3d length_scales{ 1.0, 1.0, 1.0 };

	// Grid flags
	//##############################################################
	// Describe whether a cell is active, inactive or a ghost point
	// active:		0
	// inactive:	1
	// ghost point:	2
	Eigen::MatrixXi grid_flags;

	std::map<std::pair<int, int>, Eigen::Vector2d> image_points;
	std::map<std::pair<int, int>, Eigen::Vector4d> beta_coeffs;

	Eigen::MatrixXd phi_matrix;
	Eigen::MatrixXd phi_image_point_matrix;
	Eigen::MatrixXd boundary_phi;
	Eigen::MatrixX<size_t> ghost_point_parent_sdf;

	// Weighted Least Squares Method
	//#################################################################
	// Holds copy of the data for the duration of setting the boundary
	// conditions in order to not iteratively skew the method.

	Eigen::MatrixXd m_wlsq_phi_matrix;
	std::map<std::pair<int, int>, WLSQdata> m_wlsq_data;

	// Debugging-oriented variables
	std::map<std::pair<int, int>, std::array<Eigen::Vector2d, 4>> bilinear_interp_selection;
};