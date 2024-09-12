//
// Created by 13396 on 2024/5/25.
//
#include <iostream>
#include <ducastelle.h>
#include <neighbors.h>
#include <string>
#include <xyz.h>
#include <domain.h>
#include <verlet_2.h>
#include <berendsen_thermostat.h>

void velocity_verlet(Atoms &atoms, double mass, NeighborList &neighbor_list, double timestep,
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
    verlet2_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);

    // Update neighbor list and forces
    domain.exchange_atoms(atoms);
    domain.update_ghosts(atoms, 2.0 * cutoff); // The boundary width of the EAM potential is 2 times cutoff
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
void run_nve_simulation(Atoms &atoms, double mass, double kB, double A, double xi, double p,double q, double re, NeighborList &neighbor_list,
                    double cutoff, double delta_Q, double timestep, const std::string &output_file, const std::string &output_file2, Domain &domain) {
    // Initialize file to record data
    std::ofstream data(output_file);
    std::ofstream traj(output_file2,std::ios::out);

    domain.enable(atoms);
    neighbor_list.update(atoms, cutoff);

    data << "step,total_energy,kinetic_energy,potential_energy,temperature\n";

    int step = 1000;
    int output_interval = 100;
    for (int istep1 = 0; istep1 <= step; ++istep1) {
        // Compute potential energy and forces
        // double potential_energy = ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);
        //
        // // Compute kinetic energy
        // double kinetic_energy = compute_kinetic_energy(atoms, mass);
        // double total_energy = potential_energy + kinetic_energy;

        // Velocity Verlet integration
        velocity_verlet(atoms, mass, neighbor_list, timestep, A, xi, p, q, re, cutoff, domain);

        //Compute the total potential energy without ghost atoms
        double local_potential_energy = 0.0, local_kinetic_energy = 0.0;
        for (int i = 0; i < domain.nb_local(); i++) {
            local_potential_energy += atoms.per_atom_potential_energy(i);
            local_kinetic_energy += (mass * 0.5) * atoms.velocities.col(i).square().sum();
        }

        // Sum from all processes
        double global_kinetic_energy{MPI::allreduce(local_kinetic_energy, MPI_SUM, MPI_COMM_WORLD)};
        double global_potential_energy{MPI::allreduce(local_potential_energy, MPI_SUM, MPI_COMM_WORLD)};

        double global_total_energy = global_kinetic_energy + global_potential_energy;

        //Compute temperature
        double temp_local = (2.0 * local_kinetic_energy) / (3.0 * domain.nb_local() * kB);
        double temp{MPI::allreduce(temp_local, MPI_SUM, MPI_COMM_WORLD)};

        // Only record data in main process
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            data << istep1 << "," << global_total_energy << "," << global_kinetic_energy << ","
                 << global_potential_energy << "," << temp << "\n";
        }

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
    domain.disable(atoms);
    data.close();
    traj.close();
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    try {
        Domain domain(MPI_COMM_WORLD, {50, 50, 50}, {1, 1, 1}, {0, 0, 1});
        std::string file_path_1 = "cluster_923.xyz";
        auto [names_1, positions_1]{read_xyz(file_path_1)};
        Atoms atoms_1(positions_1);

        std::string file_path_2 = "cluster_3871.xyz";
        auto [names_2, positions_2]{read_xyz(file_path_2)};
        Atoms atoms_2(positions_2);

        double cutoff = 5.0;

        // Simulation parameters
        double A = 0.2061;
        double xi = 1.790;
        double p = 10.229;
        double q = 4.036;
        double re = 2.897;
        double mass = 20405.736652;
        double kB = 8.617333262e-5;

        // Time integration parameters
        double timestep = 1; // in appropriate time units

        //Parameters of temperature
        double delta_Q = 2;       // Energy incremen

        for (int i=1; i<=2; i++){
            Atoms atoms = i==1?atoms_1:atoms_2;
            NeighborList neighbor_list;

            std::string output_file1 = "simulation_data_" + std::to_string(atoms.nb_atoms()) + ".csv";
            std::string output_file2 = "traj_" + std::to_string(atoms.nb_atoms()) + ".xyz";
            run_nve_simulation(atoms, mass, kB, A, xi, p, q, re, neighbor_list, cutoff, delta_Q,
                                timestep, output_file1, output_file2, domain);
        }
    } catch (const std::runtime_error &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    MPI_Finalize();
    return 0;
}

