
#include "ibm_application/CartGrid.h"

CartGrid::CartGrid(size_t nn) : grid_flags(nn, nn, 0), phi_matrix(nn, nn, 0)
{
    for (size_t i = 0UL; i < grid_flags.rows(); ++i) {
        for (blaze::DynamicMatrix<int, blaze::rowMajor>::Iterator it = grid_flags.begin(i); it != grid_flags.end(i); ++it) {
            auto test = *it;  // OK: Write access to the value of the element.
        }
    }
}

