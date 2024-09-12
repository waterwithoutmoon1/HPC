//
// Created by 13396 on 24-9-11.
//

#ifndef BERENDSEN_THERMOSTAT_H
#define BERENDSEN_THERMOSTAT_H
#include <atoms.h>
void apply_berendsen_thermostat(Atoms &atoms, double kB, double mass, double target_temp, double tau_T, double timestep);

double compute_kinetic_energy(const Eigen::Array3Xd &velocities, double mass);
#endif //BERENDSEN_THERMOSTAT_H
