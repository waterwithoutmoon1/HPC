//
// Created by 13396 on 2024/6/19.
//
#include <lj_neighborlist.h>
#include <atoms.h>

double lj_neighborlist(Atoms &atoms, NeighborList &neighbor_list, double epsilon, double sigma, double cutoff ) {
    double cutoff2 = cutoff * cutoff;
    double sigma2 = sigma * sigma;
    double sigma6 = sigma2 * sigma2 * sigma2;
    double cutoff_sigma6 = sigma6 / (cutoff2 * cutoff2 * cutoff2);
    double shift = 4 * epsilon * (cutoff_sigma6 * cutoff_sigma6 - cutoff_sigma6);

    double potential_energy = 0.0;
    atoms.forces.setZero();

    for (auto [i, j]: neighbor_list) {
        if (i < j) {
            Eigen::Vector3d r_ij = atoms.positions.col(i) - atoms.positions.col(j);
            double r2 = r_ij.squaredNorm();
            if (r2 < cutoff2) {
                double r2_inv = 1.0 / r2;
                double r6_inv = r2_inv * r2_inv * r2_inv;
                double sigma6_r6 = sigma6 * r6_inv;
                double lj_scalar = 24 * epsilon * (2 * sigma6_r6 * sigma6_r6 - sigma6_r6) * r2_inv;
                Eigen::Vector3d force = lj_scalar * r_ij;

                // atoms.forces.col(i) += force;
                // atoms.forces.col(j) -= force;
                atoms.forces.col(i).matrix() += force; // Apply force to atom i
                atoms.forces.col(j).matrix() -= force; // Apply opposite force to atom j

                potential_energy += 4 * epsilon * (sigma6_r6 * sigma6_r6 - sigma6_r6) - shift;
            }
        }
    }

    return potential_energy;
}