//
// Created by 13396 on 2024/5/25.
//

#include <iostream>
#include <ducastelle.h>
#include <neighbors.h>
#include <string>
#include <xyz.h>
#include <ih.cpp>

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
    double A, double xi, double p,double q, double re, double cutoff) {
    double potential_energy = ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);

    // Update positions
    atoms.velocities += 0.5 * timestep * atoms.forces / 1.0; // Assuming mass = 1.0
    atoms.positions += timestep * atoms.velocities;

    // Update neighbor list and forces
    neighbor_list.update(atoms, cutoff);
    potential_energy = ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);

    // Update velocities
    atoms.velocities += 0.5 * timestep * atoms.forces / 1.0;
}

void apply_berendsen_thermostat(Atoms &atoms, double current_temp, double target_temp, double tau_T, double timestep) {
    // double current_temp = (2.0 / 3.0) * (compute_kinetic_energy(atoms.velocities) / atoms.nb_atoms()); //(atoms.velocities.square().sum()) / (3.0 * atoms.nb_atoms());
    double lambda = sqrt(1.0 + (timestep / tau_T) * (target_temp / current_temp - 1.0));
    atoms.velocities *= lambda;
}

//Create clusters of different sizes
void create_cubic_lattice(Atoms &atoms, int num_atoms_per_side, double lattice_constant) {
    int num_atoms = num_atoms_per_side * num_atoms_per_side * num_atoms_per_side;
    Positions_t positions(3, num_atoms);
    Velocities_t velocities(3, num_atoms);
    positions.setZero();
    velocities.setZero();

    int index = 0;
    for (int x = 0; x < num_atoms_per_side; ++x) {
        for (int y = 0; y < num_atoms_per_side; ++y) {
            for (int z = 0; z < num_atoms_per_side; ++z) {
                positions(0, index) = x * lattice_constant;
                positions(1, index) = y * lattice_constant;
                positions(2, index) = z * lattice_constant;
                ++index;
            }
        }
    }
    atoms = Atoms(positions, velocities);
}

// Function to run the simulation
void run_nve_simulation(Atoms &atoms, double A, double xi, double p,double q, double re, NeighborList &neighbor_list,
                    double cutoff, double initial_temp, double final_temp, double delta_Q, double tau_T, double timestep, double relaxation_time,
                    const std::string &output_file, const std::string &output_file2) {
    // Initialize file to record data
    std::ofstream data(output_file);
    std::ofstream traj(output_file2,std::ios::out);

    data << "step,total_energy,kinetic_energy,potential_energy,temperature\n";

    double current_temp = initial_temp;
    // apply_berendsen_thermostat(atoms, current_temp, tau_T, timestep);
    std::cout << "大循环前current_temp1: " << current_temp << std::endl;

    int step = 0;

    while (current_temp < final_temp) {
        // Increase energy by rescaling velocities
        double target_temp = current_temp + delta_Q;
        double scale_factor = sqrt(target_temp / current_temp);
        atoms.velocities *= scale_factor;
        // apply_berendsen_thermostat(atoms, current_temp, target_temp, tau_T, timestep);
        int output_interval = 100;

        // Relaxation period
        int steps = static_cast<int>(relaxation_time / timestep);
        double temp_sum = 0.0;

        for (int istep = 0; istep <= steps; ++istep) {
            // Output trajectory
            if (istep % output_interval == 0) {
                write_xyz(traj, atoms);
            }

            // Velocity Verlet integration
            velocity_verlet(atoms, neighbor_list, timestep, A, xi, p, q, re, cutoff);
            temp_sum += (2.0 / 3.0) * (compute_kinetic_energy(atoms.velocities) / atoms.nb_atoms());

        }
        std::cout << "结束一次大循环后，temp之和: " << temp_sum << std::endl;
        double avg_temp = temp_sum / steps;
        // Compute potential energy and forces
        double potential_energy = ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);
        // Compute kinetic energy and temperature
        double kinetic_energy = compute_kinetic_energy(atoms.velocities);
        double total_energy = potential_energy + kinetic_energy;

        // Record data
        data << step << "," << total_energy << "," << kinetic_energy << ","
        << potential_energy << "," << avg_temp << "\n";

        // Update current temperature
        current_temp = avg_temp;
        step++;
        std::cout << "结束一次大循环后，current_temp2: " << current_temp << std::endl;
    }
    data.close();
    traj.close();
}

int main() {
    try {
        std::string file_path_1 = "./cluster_923.xyz";
        auto [names_1, positions_1]{read_xyz(file_path_1)};
        Atoms atoms_1(positions_1);

        std::string file_path_2 = "./cluster_3871.xyz";
        auto [names_2, positions_2]{read_xyz(file_path_2)};
        Atoms atoms_2(positions_2);
        // Positions_t positions;
        // Velocities_t velocities;
        // Atoms atoms(positions, velocities);

        double cutoff = 5.0;

        // Simulation parameters
        double A = 0.2061;
        double xi = 1.790;
        double p = 10.229;
        double q = 4.036;
        double re = 2.897;

        // Time integration parameters
        double relaxation_time = 1.018e-1;
        double timestep = 1.018e-3; // in appropriate time units
        double tau_T = 1.018e-1;  // Thermostat coupling constant


        //Parameters of temperature
        double initial_temp = 0.1;  // Initial temperature for the thermostat
        double final_temp = 50;    // Final target temperature
        double delta_Q = 2.5;       // Energy incremen

        for (int i=1; i<=2; i++){
            Atoms atoms = i==1?atoms_1:atoms_2;
            NeighborList neighbor_list;
            neighbor_list.update(atoms, cutoff);

            std::string output_file1 = "simulation_data_" + std::to_string(atoms.nb_atoms()) + ".csv";
            std::string output_file2 = "traj_" + std::to_string(atoms.nb_atoms()) + ".xyz";
            run_nve_simulation(atoms, A, xi, p, q, re, neighbor_list, cutoff, initial_temp, final_temp, delta_Q,
                                tau_T, timestep, relaxation_time, output_file1, output_file2);
        }

        //Simulation with different numbers of atoms
        // for (int num_atoms_per_side = 2; num_atoms_per_side <= 10; ++num_atoms_per_side) {
        //     double lattice_constant = 1.0;
        //     create_cubic_lattice(atoms, num_atoms_per_side, lattice_constant);
        //     // Initial neighbor list update
        //
        //     NeighborList neighbor_list;
        //     neighbor_list.update(atoms, cutoff);
        //
        //     std::string output_file = "simulation_data_" + std::to_string(atoms.nb_atoms()) + ".csv";
        //     run_nve_simulation(atoms, A, xi, p, q, re, neighbor_list, cutoff, initial_temp, final_temp, tau_T, delta_Q,
        //                     timestep, relaxation_time, output_file);
        // }
    } catch (const std::runtime_error &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}




