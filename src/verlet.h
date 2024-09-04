//
// Created by 13396 on 2024/5/25.
//

#ifndef VERLET_H
#define VERLET_H

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep, double mass);
void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep, double mass);

#endif //VERLET_H
