
#include "ibm_application/CartGrid.h"

CartGrid::CartGrid(size_t nn, double phi_0) : grid_flags(nn, nn, 0), image_points(nn,nn), phi_matrix(nn, nn, phi_0)
{
    h = length_scales.at(1) / static_cast<double>(nn-1);
}

void CartGrid::UpdateGrid()
{
    for (auto const& [hash, sdf] : immersed_boundaries)
    {
        // For each boundary (given that the boundary data exists) we update the grid flags
        if (sdf)
        {
            for (size_t i = 0UL; i < grid_flags.rows(); ++i) {
                for (size_t j = 0UL; j < grid_flags.columns(); ++j) {

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
                            blaze::StaticVector<double, 2L> node{ x, y };

                            auto image_point = node + normal * std::abs(sdf->SignedDistanceFunction(x, y)) * 2.0;
                            image_points(i, j) = image_point;
                        }
                        else {
                            grid_flags(i, j) = 1;
                        }
                    }
                    else {
                        grid_flags(i, j) = 0;
                    }
                }
            }
        }
    }

    std::cout << grid_flags;
    std::cout << image_points;
}
