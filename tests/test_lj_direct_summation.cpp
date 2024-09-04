
#include <gtest/gtest.h>
#include <atoms.h>
#include "lj_direction_summation.h"

TEST(LJDirectSummationTest, Forces) {
    constexpr int nb_atoms = 10;
    constexpr double epsilon = 0.7;  // choose different to 1 to pick up missing factors
    constexpr double sigma = 0.3;
    constexpr double delta = 0.0001;  // difference used for numerical (finite difference) computation of forces
    Positions_t positions(3,nb_atoms);

    Atoms atoms(positions);
    atoms.positions.setRandom();  // random numbers between -1 and 1

    // compute and store energy of the indisturbed configuration

    double e0{lj_direct_summation(atoms, epsilon, sigma)};
    Forces_t forces0{atoms.forces};

    // loop over all atoms and compute forces from a finite differences approximation
    for (int i{0}; i < nb_atoms; ++i) {
        // loop over all Cartesian directions
        for (int j{0}; j < 3; ++j) {
            // move atom to the right
            atoms.positions(j, i) += delta;
            double eplus{lj_direct_summation(atoms, epsilon, sigma)};
            // move atom to the left
            atoms.positions(j, i) -= 2 * delta;
            double eminus{lj_direct_summation(atoms, epsilon, sigma)};
            // move atom back to original position
            atoms.positions(j, i) += delta;

            // finite-differences forces
            double fd_force{-(eplus - eminus) / (2 * delta)};

            std::cout << "fd_force:" << fd_force << std::endl;
            std::cout << "forces0:" << forces0(j,i) << std::endl;

            // check whether finite-difference and analytic forces agree
            if (abs(forces0(j, i)) > 1e-10) {
                EXPECT_NEAR(abs(fd_force - forces0(j, i)) / forces0(j, i), 0, 1e-5);
            } else {
                EXPECT_NEAR(fd_force, forces0(j, i), 1e-10);
            }
        }
    }
}