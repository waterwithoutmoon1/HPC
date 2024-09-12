//
// Created by 13396 on 2024/5/25.
//

#include <iostream>
#include <lj_direction_summation.cpp>

#include <string>
#include <xyz.h>
#include <verlet_2.h>

double compute_kinetic_energy(const Eigen::Array3Xd &velocities, double mass) {
    double kinetic_energy = 0.0;
    size_t num_atoms = velocities.cols();
    for (int i = 0; i < num_atoms; ++i) {
        kinetic_energy += 0.5 * mass * velocities.col(i).matrix().squaredNorm();
    }
    return kinetic_energy;
}

void velocity_verlet(Atoms &atoms, double timestep, double epsilon, double sigma, double mass) {
    size_t num_atoms = atoms.nb_atoms();
    double potential_energy = lj_direct_summation(atoms, epsilon, sigma);
    verlet2_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);

    potential_energy = lj_direct_summation(atoms, epsilon, sigma);
    verlet2_step2(atoms.velocities, atoms.forces, timestep, mass);
}

int main() {
try{
    std::string file_path = "lj54.xyz";
    auto [names, positions, velocities]{read_xyz_with_velocities(file_path)};
    Atoms atoms(positions, velocities);

    double epsilon = 1.0;
    double sigma = 1.0;
    double mass = 1.0;
    double dt = 0.01;
    double total_time = 100.0;
    int num_steps = static_cast<int>(total_time / dt);
    int output_interval = 100;

    std::ofstream traj("./traj.xyz",std::ios::out);
    std::ofstream data("./simulation_data.csv",std::ios::out);

    data << "step,total_energy,kinetic_energy,potential_energy\n";
    for (int step = 0; step < num_steps; ++step) {
        if (step % output_interval == 0) {
            write_xyz(traj, atoms);
        }

        double kinetic_energy = compute_kinetic_energy(atoms.velocities, mass);
        double potential_energy = lj_direct_summation(atoms, epsilon, sigma);
        double total_energy = kinetic_energy + potential_energy;

        std::cout << "Step: " << step
                  << " Total Energy: " << total_energy
                  << " Kinetic Energy: " << kinetic_energy
                  << " Potential Energy: " << potential_energy << std::endl;

        data << step << "," << total_energy << "," << kinetic_energy << ","
                 << potential_energy << "\n";

        velocity_verlet(atoms, dt, epsilon, sigma, mass);
    }
    traj.close();
    data.close();
} catch (const std::runtime_error &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}