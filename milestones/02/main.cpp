//
// Created by 13396 on 2024/5/25.
//
#include <Eigen/Dense>
#include <iostream>
#include <verlet.h>

int main() {

  int nb_steps = 10;
  double x, y, z, vx, vy, vz = 0;
  double fx, fy ,fz = 1;
  double timestep = 0.5;
  double mass = 0.2;

  for (int i = 0; i < nb_steps; ++i) {
    double vx1 = vx;
    double vy1 = vy;
    double vz1 = vz;
    std::cout << "Step: " << i << std::endl;

    verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, timestep, mass);

    fx = mass * (vx - vx1)/ timestep;
    fy = mass * (vy - vy1)/ timestep;
    fz = mass * (vz - vz1)/ timestep;

    verlet_step2(vx, vy, vz, fx, fy, fz, timestep, mass);
      std::cout << "x" << x << std::endl;
  }
  return 0;
}