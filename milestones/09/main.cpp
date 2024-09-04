//
// Created by 13396 on 2024/5/25.
//
#include <iostream>
#include <ducastelle.h>
#include <neighbors.h>
#include <string>
#include <xyz.h>
#include <domain.h>

// Function to compute kinetic energy
double compute_kinetic_energy(const Eigen::Array3Xd &velocities) {
    double kinetic_energy = 0.0;
    size_t num_atoms = velocities.cols();
    for (int i = 0; i < num_atoms; ++i) {
        kinetic_energy += 0.5 * 1 * velocities.col(i).matrix().squaredNorm();
    }
    return kinetic_energy;
}

void velocity_verlet(Atoms &atoms, NeighborList &neighbor_list, double timestep,
    double A, double xi, double p,double q, double re, double cutoff, Domain &domain) {

    double potential_energy = 0.0;
    try {
        std::cout << "Number of atoms: " << atoms.nb_atoms() << std::endl;
        std::cout << "Neighbor list size: " << neighbor_list.nb_neighbors() << std::endl;
        potential_energy = ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);
    } catch (const std::exception &e) {
        std::cerr << "Error in ducastelle: " << e.what() << std::endl;
    }
    // Update positions
    atoms.velocities += 0.5 * timestep * atoms.forces / 1.0; // Assuming mass = 1.0
    atoms.positions += timestep * atoms.velocities;

    // Update neighbor list and forces
    neighbor_list.update(atoms, cutoff);
    domain.exchange_atoms(atoms);
    domain.update_ghosts(atoms, 2.0 * cutoff); // EAM势能的边界宽度为2倍cutoff
    potential_energy = ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);

    // Update velocities
    atoms.velocities += 0.5 * timestep * atoms.forces / 1.0;
}

// Function to compute stress tensor
Eigen::Array3Xd compute_cauchy_stress(const Atoms &atoms, const NeighborList &neighbor_list, double volume) {
    Eigen::Matrix3d stress_1 = Eigen::Matrix3d::Zero();
    for (auto [i, j] : neighbor_list) {
        if (i < j) {
            Eigen::Vector3d r_ij = atoms.positions.col(i) - atoms.positions.col(j);
            Eigen::Vector3d f_ij = atoms.forces.col(i) - atoms.forces.col(j);
            stress_1 += r_ij * f_ij.transpose();
        }
    }
    Eigen::Array3Xd stress = stress_1.array();
    return stress / volume;
}


// Function to run the simulation
void run_rve_simulation(Atoms &atoms, double A, double xi, double p, double q, double re, NeighborList &neighbor_list,
                        double cutoff, double initial_temp, double delta_Q, double timestep, double relaxation_time, const std::string &output_file, const std::string &output_file2,
                        Domain &domain, double strain_rate, double max_strain) {

    // Initialize file to record data
    std::ofstream data(output_file);
    std::ofstream traj(output_file2, std::ios::out);

    domain.enable(atoms);

    data << "step,total_energy,kinetic_energy,potential_energy,temperature,strain,stress_xx,stress_yy,stress_zz\n";

    double current_temp = initial_temp;
    double initial_length = domain.domain_length(2); // Initial length along z-axis
    Eigen::Array3d domain_length = domain.domain_length();
    std::cout << "initial_length " << domain_length[2] << std::endl;
    int step = 0;
    int output_interval = 100;

    // Strain loop
    double strain = 0.0;
    while (strain <= max_strain) {
        // Increase energy by rescaling velocities
        double target_temp = current_temp + delta_Q;
        double scale_factor = sqrt(target_temp / current_temp);
        atoms.velocities *= scale_factor;

        // Relaxation period
        int steps = static_cast<int>(relaxation_time / timestep);
        double temp_sum = 0.0;
        for (int istep = 1; istep <= steps; ++istep) {
            // Output trajectory
            if (istep % output_interval == 0) {
                write_xyz(traj, atoms);
            }

            // Velocity Verlet integration
            velocity_verlet(atoms, neighbor_list, timestep, A, xi, p, q, re, cutoff, domain);
            temp_sum += (2.0 / 3.0) * (compute_kinetic_energy(atoms.velocities) / atoms.nb_atoms());
        }

        // Compute average temperature
        double avg_temp = temp_sum / steps;

        // Compute potential energy and forces
        double potential_energy = ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);

        // Compute kinetic energy and temperature
        double kinetic_energy = compute_kinetic_energy(atoms.velocities);
        double total_energy = potential_energy + kinetic_energy;

        // Compute Cauchy stress tensor
        Eigen::Array3Xd stress_tensor = compute_cauchy_stress(atoms, neighbor_list, domain.volume());

        // Record data
        data << step << "," << total_energy << "," << kinetic_energy << "," << potential_energy << ","
             << avg_temp << "," << strain << ","
             << stress_tensor(0, 0) << "," << stress_tensor(1, 1) << "," << stress_tensor(2, 2) << "\n";

        // Stretch the domain
        strain += strain_rate * timestep;
        domain_length[2] = initial_length * (1.0 + strain); // Update length along the z-axis
        domain.scale(atoms, domain_length); // Rescale the domain and atomic positions
        neighbor_list.update(atoms, cutoff); // Update neighbor list after scaling
        std::cout << "strain: " << domain.domain_length(2) << std::endl;
        // Update current temperature
        current_temp = avg_temp;
        step++;
    }

    domain.disable(atoms);
    data.close();
    traj.close();
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    try {
        Domain domain(MPI_COMM_WORLD, {20, 20, 100}, {1, 1, 1}, {0, 0, 1});
        std::string file_path_1 = "./whisker_small.xyz";
        auto [names_1, positions_1]{read_xyz(file_path_1)};
        Atoms atoms(positions_1);

        double cutoff = 2.5;

        // Simulation parameters
        double A = 0.2061;
        double xi = 1.790;
        double p = 10.229;
        double q = 4.036;
        double re = 2.897;

        // Time integration parameters
        double relaxation_time = 1.018e-1;
        double timestep = 1.018e-3; // in appropriate time units
        // double tau_T = 1.018e-1;  // Thermostat coupling constant

        //Parameters of temperature
        double initial_temp = 0.1;  // Initial temperature for the thermostat
        // double final_temp = 50;    // Final target temperature
        double delta_Q = 2;       // Energy incremen
        double strain_rate = 10;
        double max_strain = 0.2;

        NeighborList neighbor_list;
        neighbor_list.update(atoms, cutoff);

        std::string output_file1 = "simulation_data_" + std::to_string(atoms.nb_atoms()) + ".csv";
        std::string output_file2 = "traj_" + std::to_string(atoms.nb_atoms()) + ".xyz";
        run_rve_simulation(atoms, A, xi, p, q, re, neighbor_list, cutoff, initial_temp, delta_Q,
                                timestep, relaxation_time, output_file1, output_file2, domain, strain_rate, max_strain);

    } catch (const std::runtime_error &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    MPI_Finalize();
    return 0;
}
