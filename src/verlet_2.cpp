//
// Created by 13396 on 2024/6/4.
//
#include "verlet_2.h"

void verlet2_step1(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                  const Eigen::Array3Xd &forces, double timestep, double mass) {
    velocities += 0.5 * forces * timestep / mass;
    positions += velocities * timestep;
}

void verlet2_step2(Eigen::Array3Xd &velocities, const Eigen::Array3Xd &forces, double timestep, double mass) {
    velocities += 0.5 * forces * timestep / mass;
}
