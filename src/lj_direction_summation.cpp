//
// Created by 13396 on 2024/6/5.
//
#include <atoms.h>
#include <lj_direction_summation.h>

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma){
    atoms.forces.setZero();
    double potential_energy = 0.0;
    size_t num_atoms = atoms.nb_atoms();
    for (int i = 0; i < num_atoms; ++i) {
        for (int j = i + 1; j < num_atoms; ++j) {
            // Distance vector between atoms i and j
            Eigen::Vector3d r_ij;
            // Squared distance
            double r2 = 0.0;

            // Compute distance vector and squared distance
            for (int k = 0; k < 3; ++k) {
                r_ij[k] = atoms.positions(k,i) - atoms.positions(k,j);
                r2 += r_ij[k] * r_ij[k];
            }

            // Compute distance vector and squared distance
            double r6 = r2 * r2 * r2;
            double r12 = r6 * r6;
            double sigma6 = std::pow(sigma, 6);
            double sigma12 = sigma6 * sigma6;

            double inv_r6 = sigma6 / r6;
            double inv_r12 = sigma12 / r12;

            // Lennard-Jones potential for the pair
            double potential_pair = 4.0 * epsilon * (inv_r12 - inv_r6);
            potential_energy += potential_pair;

            // Magnitude of the force
            double force_magnitude = 24.0 * epsilon * (2 * inv_r12 - inv_r6) / r2;

            // Update forces
            for (int k = 0; k < 3; ++k) {
                double force_component = force_magnitude * r_ij[k];
                atoms.forces(k,i) += force_component;
                atoms.forces(k,j) -= force_component;
            }
        }
    }
    return potential_energy;

}