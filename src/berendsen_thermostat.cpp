//
// Created by 13396 on 24-9-11.
//
#include <berendsen_thermostat.h>

double compute_kinetic_energy(const Eigen::Array3Xd &velocities, double mass) {
    double kinetic_energy = 0.0;
    size_t num_atoms = velocities.cols();
    for (int i = 0; i < num_atoms; i++) {
        kinetic_energy += 0.5 * mass * velocities.col(i).matrix().squaredNorm();
    }
    return kinetic_energy;
}

void apply_berendsen_thermostat(Atoms &atoms, double kB, double mass, double target_temp, double tau_T, double timestep) {
    double current_temp = (2.0 / 3.0) * (compute_kinetic_energy(atoms.velocities, mass) / (atoms.nb_atoms() * kB)); //(atoms.velocities.square().sum()) / (3.0 * atoms.nb_atoms());
    double lambda = sqrt(1.0 + (timestep / tau_T) * (target_temp / current_temp - 1.0));
    atoms.velocities *= lambda;
}