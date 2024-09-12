//
// Created by 13396 on 2024/5/25.
//

#include <iostream>
#include <ducastelle.h>
#include <neighbors.h>
#include <string>
#include <xyz.h>
#include <verlet_2.h>
#include <berendsen_thermostat.h>

void velocity_verlet(Atoms &atoms, double mass, NeighborList &neighbor_list, double timestep,
    double A, double xi, double p,double q, double re, double cutoff) {
    double potential_energy = ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);
    // Update positions
    verlet2_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);

    // Update neighbor list and forces
    neighbor_list.update(atoms, cutoff);
    potential_energy = ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);

    // Update velocities
    verlet2_step2(atoms.velocities, atoms.forces, timestep, mass);
}

void rescale_velocities(Atoms& atoms, double delta_Q, double mass){
    double current_kinetic_energy = compute_kinetic_energy(atoms.velocities, mass);
    double new_kinetic_energy = current_kinetic_energy + delta_Q;
    double scale_factor = sqrt(new_kinetic_energy / current_kinetic_energy);
    atoms.velocities *= scale_factor;
}

// Function to run the simulation
void run_nve_simulation(Atoms &atoms, double mass, double A, double xi, double p,double q, double re, NeighborList &neighbor_list,
                    double cutoff, double delta_Q, double kB, double timestep, const std::string &output_file, const std::string &output_file2) {
    // Initialize file to record data
    std::ofstream data(output_file);
    std::ofstream traj(output_file2,std::ios::out);

    data << "step,total_energy,kinetic_energy,potential_energy,temperature\n";

    int step = 1000;
    int output_interval = 100;
    for (int istep1 = 0; istep1 <= step; istep1++) {
        // Compute potential energy and forces
        double potential_energy = ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);

        // Compute kinetic energy
        double kinetic_energy = compute_kinetic_energy(atoms.velocities, mass);
        double total_energy = potential_energy + kinetic_energy;

        //Compute temperature
        double temp = (2.0 * kinetic_energy) / (3.0 *  atoms.nb_atoms() * kB);

        // Velocity Verlet integration
        velocity_verlet(atoms, mass, neighbor_list, timestep, A, xi, p, q, re, cutoff);

        // Record data
        data << istep1 << "," << total_energy << "," << kinetic_energy << ","
        << potential_energy << "," << temp << "\n";

        // Output trajectory
        if (istep1 % output_interval == 0) {
            write_xyz(traj, atoms);
        }

        // Increase energy by rescaling velocities
        if (istep1 % 10 == 0 && istep1 != 0){
            rescale_velocities(atoms, delta_Q, mass);
        }

        // Update current temperature
        std::cout << "After a loop，current_temp2: " << temp << std::endl;
    }
    data.close();
    traj.close();
}

int main() {
    try {
        std::string file_path_1 = "./cluster_55.xyz";
        auto [names_1, positions_1]{read_xyz(file_path_1)};
        Atoms atoms_1(positions_1);

        std::string file_path_2 = "./cluster_147.xyz";
        auto [names_2, positions_2]{read_xyz(file_path_2)};
        Atoms atoms_2(positions_2);

        double cutoff = 5.0;

        // Simulation parameters
        double A = 0.2061;
        double xi = 1.790;
        double p = 10.229;
        double q = 4.036;
        double re = 2.897;

        // Time integration parameters
        double timestep = 1.0; // in appropriate time units

        //Parameters of temperature
        double delta_Q = 1;       // Energy incremen

        double mass = 20405.736652;
        double kB = 8.617333262e-5;

        for (int i=1; i<=2; i++){
            Atoms atoms = i==1?atoms_1:atoms_2;
            NeighborList neighbor_list;
            neighbor_list.update(atoms, cutoff);

            std::string output_file1 = "simulation_data_" + std::to_string(atoms.nb_atoms()) + ".csv";
            std::string output_file2 = "traj_" + std::to_string(atoms.nb_atoms()) + ".xyz";
            run_nve_simulation(atoms, mass, A, xi, p, q, re, neighbor_list, cutoff, delta_Q,
                                kB, timestep, output_file1, output_file2);
        }
    } catch (const std::runtime_error &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}




