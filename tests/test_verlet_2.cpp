//
// Created by 13396 on 2024/6/4.
//
#include "verlet_2.h"
#include <gtest/gtest.h>
#include <iostream>
#include <Eigen/Dense>

// Demonstrate some basic assertions.
TEST(VerletTest2, Loops) {
    int nb_steps = 10;
    int nb_atoms = 10;
    using Positions_t = Eigen::Array3Xd;
    using Velocities_t = Eigen::Array3Xd;
    using Forces_t = Eigen::Array3Xd;
    Positions_t positions(3, nb_atoms);
    Velocities_t velocities(3, nb_atoms);
    Forces_t forces(3, nb_atoms);
    double timestep = 0.5;
    double mass = 0.2;

    for (int i = 0; i < nb_steps; ++i) {
        Eigen::Array3Xd velocities_1 = velocities;
        std::cout << "Step: " << i << std::endl;

        verlet2_step1(positions, velocities, forces, timestep, mass);

        forces = mass * (velocities - velocities_1)/ timestep;

        verlet2_step2(velocities, forces, timestep, mass);

        std::cout << "positions_x" << positions.row(0) << std::endl;
        std::cout << "positions_y" << positions.row(1) << std::endl;
        std::cout << "positions_z" << positions.row(2) << std::endl;
        std::cout << "velocities_x" << velocities.row(0) << std::endl;
        std::cout << "velocities_y" << velocities.row(1) << std::endl;
        std::cout << "velocities_z" << velocities.row(2) << std::endl;
    }


}