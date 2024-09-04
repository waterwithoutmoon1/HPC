//
// Created by 13396 on 2024/5/25.
//
#include "verlet.h"


void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep, double mass) {
    vx += 0.5 * fx * timestep / mass;
    x += vx * timestep;
    vy += 0.5 * fy * timestep / mass;
    y += vy * timestep;
    vz += 0.5 * fz * timestep / mass;
    z += vz * timestep;
}

void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep, double mass) {
    vx += 0.5 * fx * timestep / mass;
    vy += 0.5 * fy * timestep / mass;
    vz += 0.5 * fz * timestep / mass;
}
