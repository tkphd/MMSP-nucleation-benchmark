// mmsp2norm.cpp
// Compute L2-norm of two grids from the same time, with 2:1 resolution
// Questions/comments to trevor.keller@nist.gov (Trevor Keller)

#include <cmath>
#include <iostream>

#include "MMSP.hpp"
#include "nucleation.hpp"

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cout << "Error: " << argv[0] << " requires two MMSP grids as input." << std::endl;
        std::cout << "Usage: " << argv[0] << " grid_dx.dat grid_2dx.dat" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    GRID2D gridA(argv[1]);
    GRID2D gridB(argv[2]);

    if (MMSP::dx(gridB) < MMSP::dx(gridA)) {
        std::cout << "Error: Please supply the more-refined grid first." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    GRID2D gridE(gridB);
    for (int n = 0; n < MMSP::nodes(gridE); n++)
        gridE(n) = 0.0;

    MMSP::vector<int> xA = MMSP::position(gridA, 0);
    MMSP::vector<int> xB = MMSP::position(gridB, 0);

    const int stride = std::ceil(MMSP::dx(gridB) / MMSP::dx(gridA));

    double e2 = 0.0;

    while (xB[0] < MMSP::x1(gridB) && xA[0] < MMSP::x1(gridA)) {
        xB[1] = y0(gridB);
        xA[1] = y0(gridA);
        while (xB[1] < MMSP::y1(gridB) && xA[1] < MMSP::y1(gridA)) {
            const double delta_phi = gridA(xA) - gridB(xB);
            const double err = std::sqrt(delta_phi * delta_phi);
            gridE(xB) = err;
            e2 += err;

            xB[1] += 1;
            xA[1] += stride;
        }
        xB[0] += 1;
        xA[0] += stride;
    }

    double dV = 1.0;
    for (int d = 0; d < 2; d++)
        dV *= MMSP::dx(gridB, d);

    e2 *= dV;

    std::cout << MMSP::dx(gridA) << ',' << MMSP::dx(gridB) << ',' << e2 << std::endl;

    MMSP::output(gridE, "error.dat");

    return EXIT_SUCCESS;
}
