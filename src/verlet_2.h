//
// Created by 13396 on 2024/6/4.
//

#ifndef VERLET_2_H
#define VERLET_2_H
#include <Eigen/Dense>

void verlet2_step1(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                  const Eigen::Array3Xd &forces, double timestep, double mass);

void verlet2_step2(Eigen::Array3Xd &velocities, const Eigen::Array3Xd &forces, double timestep, double mass);

#endif //VERLET_2_H
