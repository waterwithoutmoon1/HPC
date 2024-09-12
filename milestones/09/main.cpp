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
// // Function to compute kinetic energy
// double compute_kinetic_energy(Atoms &atoms, double mass, Domain &domain) {
//     double kinetic_energy = 0.0;
//     size_t num_atoms = atoms.velocities.cols();
//     for (int i = 0; i < domain.nb_local(); ++i) {
//         kinetic_energy += 0.5 * mass * atoms.velocities.col(i).matrix().squaredNorm();
//     }
//     return kinetic_energy;
// }
//
// void velocity_verlet(Atoms &atoms, double mass, NeighborList &neighbor_list, double timestep,
//     double A, double xi, double p,double q, double re, double cutoff, Domain &domain) {
//     double potential_energy = 0.0;
//     try {
//         std::cout << "Number of atoms: " << atoms.nb_atoms() << std::endl;
//         std::cout << "Neighbor list size: " << neighbor_list.nb_neighbors() << std::endl;
//         potential_energy = ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);
//     } catch (const std::exception &e) {
//         std::cerr << "Error in ducastelle: " << e.what() << std::endl;
//     }
//     // Update positions
//     verlet2_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);
//
//     // Update neighbor list and forces
//     domain.exchange_atoms(atoms);
//     domain.update_ghosts(atoms, 2.0 * cutoff); // EAM势能的边界宽度为2倍cutoff
//     neighbor_list.update(atoms, cutoff);
//     potential_energy = ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);
//
//     // Update velocities
//     verlet2_step2(atoms.velocities, atoms.forces, timestep, mass);
// }
//
// void apply_berendsen_thermostat(Atoms &atoms, double kB, double mass, double target_temp, double tau_T, double timestep, Domain &domain) {
//     double current_temp = (2.0 / 3.0) * (compute_kinetic_energy(atoms, mass, domain) / (atoms.nb_atoms() * kB));
//     double lambda = sqrt(1.0 + (timestep / tau_T) * (target_temp / current_temp - 1.0));
//     std::cout <<"lamda:"<< lambda << std::endl;
//     for (int i = 0; i < domain.nb_local(); ++i) {
//         atoms.velocities.col(i) *= lambda;
//     }
// }
//
// // Function to run the simulation
// void run_rve_simulation(Atoms &atoms, double mass, double kB, double A, double xi, double p, double q, double re, NeighborList &neighbor_list,
//                         double cutoff, double target_temp, double tau_T, double timestep, const std::string &output_file, const std::string &output_file2,
//                         Domain &domain, double max_scale) {
//     // Initialize file to record data
//     std::ofstream data(output_file);
//     std::ofstream traj(output_file2,std::ios::out);
//
//     domain.enable(atoms);
//     neighbor_list.update(atoms, cutoff);
//
//     data << "step,total_energy,kinetic_energy,potential_energy,temperature,strain,stress\n";
//
//     long scale_interval = 1000;
//     double scale_amount = 0.5;
//     double domain_size_z = domain.domain_length(2);
//     // FOR SMALL WHISKER
//     double area = 1600;
//
//     // scale rate (unit = 1/fs)
//     double scale_rate = scale_amount / (scale_interval * domain_size_z);
//     std::cout << "scale_rate:" << scale_rate << std::endl;
//
//     double total_scaled = 0;
//
//     int i = 0;
//     while (total_scaled < max_scale) {
//         total_scaled += scale_amount;
//         domain.scale(atoms, Eigen::Vector3d(domain.domain_length(0), domain.domain_length(1), domain_size_z + total_scaled));
//         neighbor_list.update(atoms, cutoff);
//
//         velocity_verlet(atoms, mass, neighbor_list, timestep, A,  xi,  p, q,  re, cutoff,  domain);
//
//         //Compute temperature
//         double temp_local = (2.0 * compute_kinetic_energy(atoms, mass, domain)) / (3.0 *  atoms.nb_atoms() * kB);
//         std::cout << "temp_local:" << temp_local << std::endl;
//         //Compute the total potential energy without ghost atoms
//         double local_potential_energy = 0.0, local_kinetic_energy = 0.0;
//         for (int i = 0; i < domain.nb_local(); i++) {
//             local_potential_energy += atoms.per_atom_potential_energy(i);
//             local_kinetic_energy += (mass * 0.5) * atoms.velocities.col(i).square().sum();
//         }
//
//         // Sum from all processes
//         double global_kinetic_energy{MPI::allreduce(local_kinetic_energy, MPI_SUM, MPI_COMM_WORLD)};
//         double global_potential_energy{MPI::allreduce(local_potential_energy, MPI_SUM, MPI_COMM_WORLD)};
//         double global_total_energy = global_kinetic_energy + global_potential_energy;
//
//         // compute the forces on ghost atoms outside the boundaries
//         double local_force1 = 0.0, local_force2 = 0.0;
//
//         for (int j = domain.nb_local(); j < atoms.nb_atoms(); ++j) {
//             // Ghost atoms to the bottom of the whisker
//             if (atoms.positions(2, j) < atoms.positions.row(2).minCoeff())
//                 local_force1 += atoms.forces(2, j);
//
//             // Ghost atoms to the top of the whisker
//             if (atoms.positions(2, j) > domain_size_z + total_scaled)
//                 local_force2 += atoms.forces(2, j);
//         }
//
//         // Sum over all processes
//         double global_force1{MPI::allreduce(local_force1, MPI_SUM, MPI_COMM_WORLD)};
//         double global_force2{MPI::allreduce(local_force2, MPI_SUM, MPI_COMM_WORLD)};
//
//         apply_berendsen_thermostat(atoms, kB, mass, target_temp, tau_T, timestep,domain);
//         std::cout << "velocity:" << atoms.velocities.col(0) << std::endl;
//
//         double temp = MPI::allreduce(temp_local, MPI_SUM, MPI_COMM_WORLD) / domain.size();
//         std::cout << "temp:" << temp << std::endl;
//         if (domain.rank() == 0) {
//             double strain = total_scaled / domain_size_z;
//             double stress = (global_force1 - global_force2) / area;
//             data << i << "," << global_total_energy << "," << global_kinetic_energy << "," << global_potential_energy << ","
//             << temp << "," << strain << "," << stress << "\n";
//
//             if (i % 5 == 0) write_xyz(traj, atoms);
//
//             std::cout << i << std::endl;
//         }
//         i++;
//     }
//     domain.disable(atoms);
//     data.close();
//     traj.close();
// }
// int main(int argc, char *argv[]) {
//     MPI_Init(&argc, &argv);
//     try {
//
//         std::string file_path_1 = "./whisker_small.xyz";
//         auto [names_1, positions_1]{read_xyz(file_path_1)};
//         Atoms atoms(positions_1);
//
//         // Set up domains
//         double whisker_x = atoms.positions.row(0).maxCoeff();
//         double whisker_y = atoms.positions.row(1).maxCoeff();
//         double whisker_z = atoms.positions.row(2).maxCoeff() + 1;
//         std::cout << "Whisker dimensions: " << whisker_x << ", " << whisker_y << ", " << whisker_z << "\n";
//
//         Domain domain(MPI_COMM_WORLD, {whisker_x, whisker_y, whisker_z}, {1, 1, 1}, {0, 0, 1});
//
//         double cutoff = 10.0;
//
//         // Simulation parameters
//         double A = 0.2061;
//         double xi = 1.790;
//         double p = 10.229;
//         double q = 4.036;
//         double re = 2.897;
//         double mass = 20405.736652;
//         double kB = 8.617333262e-5;
//
//         // Time integration parameters
//         double timestep = 1; // in appropriate time units
//         double tau_T = 1000;
//
//         //Parameters of temperature
//         double target_temp = 1;
//         double max_scale = 100.0;
//
//         NeighborList neighbor_list;
//
//         std::string output_file1 = "simulation_data_" + std::to_string(atoms.nb_atoms()) + ".csv";
//         std::string output_file2 = "traj_" + std::to_string(atoms.nb_atoms()) + ".xyz";
//         run_rve_simulation(atoms, mass, kB, A, xi, p, q, re, neighbor_list, cutoff, target_temp, tau_T,
//                                 timestep, output_file1, output_file2, domain, max_scale);
//
//     } catch (const std::runtime_error &e) {
//         std::cerr << "Error: " << e.what() << std::endl;
//     }
//     MPI_Finalize();
//     return 0;
// }


void velocity_verlet(Atoms &atoms, double mass, NeighborList &neighbor_list, double timestep,
    double A, double xi, double p,double q, double re, double cutoff, Domain &domain) {
    double potential_energy = 0.0;
    potential_energy = ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);

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


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    auto [names, positions]{read_xyz("whisker_large.xyz")};
    Atoms atoms{positions};

    std::string output_file1 = "simulation_data_" + std::to_string(atoms.nb_atoms()) + ".csv";
    std::string output_file2 = "traj_" + std::to_string(atoms.nb_atoms()) + ".xyz";
    std::ofstream data(output_file1);
    std::ofstream traj(output_file2);


    data << "step,total_energy,kinetic_energy,potential_energy,temperature,strain,stress\n";

    // Set up domain
    double whisker_x = atoms.positions.row(0).maxCoeff();
    double whisker_y = atoms.positions.row(1).maxCoeff();
    double whisker_z = atoms.positions.row(2).maxCoeff() + 1;
    double whisker_z_min = atoms.positions.row(2).minCoeff();
    std::cout << "Whisker dimensions: " << whisker_x << ", " << whisker_y << ", " << whisker_z << "\n";

    Domain domain(MPI_COMM_WORLD, {whisker_x, whisker_y, whisker_z},
                  {1, 1, MPI::comm_size(MPI_COMM_WORLD)}, {0, 0, 1});

    //Parameters of atoms
    double A = 0.2061;
    double xi = 1.790;
    double p = 10.229;
    double q = 4.036;
    double re = 2.897;
    constexpr double mass = 20405.736652;
    constexpr double cutoff = 10.0;
    constexpr double kB = 8.617333262e-5;

    constexpr double relaxation_time = 1000;
    constexpr double target_temp = 20;

    NeighborList neighbor_list;

    domain.enable(atoms);
    neighbor_list.update(atoms, cutoff);

    double timestep = 1;

    long scale_start = 1000;//Start of record
    long scale_interval = 100;//The interval of scaling
    double scale_amount = 0.5;

    // For area
    double area = 1600;

    // scale rate (unit = 1/fs)
    double scale_rate = scale_amount / (scale_interval * whisker_z);
    std::cout << "strain rate:" <<scale_rate << std::endl;

    double max_scale_amount = 20.0;
    double total_scaled = 0;

    int i = 0;
    // The loop goes on until scaled to the max value
    while (total_scaled < max_scale_amount) {
        if (i >= scale_start && i % scale_interval == 0) {
            total_scaled += scale_amount;
            domain.scale(atoms, Eigen::Vector3d(whisker_x, whisker_y, whisker_z + total_scaled));
            neighbor_list.update(atoms, cutoff);
        }
        velocity_verlet(atoms, mass, neighbor_list, timestep, A,  xi,  p, q,  re, cutoff,  domain);

        apply_berendsen_thermostat(atoms, kB, mass, target_temp, relaxation_time, timestep);

        // compute the forces on ghost atoms outside the boundaries
        double local_force1 = 0.0, local_force2 = 0.0;

        for (int j = domain.nb_local(); j < atoms.nb_atoms(); j++) {
            // Ghost atoms to the bottom of the whisker
            if (atoms.positions(2, j) < whisker_z_min)
                local_force1 += atoms.forces(2, j);

            // Ghost atoms to the top of the whisker
            if (atoms.positions(2, j) > whisker_z + total_scaled)
                local_force2 += atoms.forces(2, j);
        }

        // Sum over all processes
        double global_force1{MPI::allreduce(local_force1, MPI_SUM, MPI_COMM_WORLD)};
        double global_force2{MPI::allreduce(local_force2, MPI_SUM, MPI_COMM_WORLD)};

        //Compute the total potential energy without ghost atoms
        double local_potential_energy = 0.0, local_kinetic_energy = 0.0;
        for (int k = 0; k < domain.nb_local(); k++) {
            local_potential_energy += atoms.per_atom_potential_energy(k);
            local_kinetic_energy += (mass * 0.5) * atoms.velocities.col(k).square().sum();
        }

        // Sum from all processes
        double global_kinetic_energy{MPI::allreduce(local_kinetic_energy, MPI_SUM, MPI_COMM_WORLD)};
        double global_potential_energy{MPI::allreduce(local_potential_energy, MPI_SUM, MPI_COMM_WORLD)};

        double global_total_energy = global_kinetic_energy + global_potential_energy;

        //Compute temperature
        double temp_local = (2.0 * local_kinetic_energy) / (3.0 *  domain.nb_local() * kB);
        double temp{MPI::allreduce(temp_local, MPI_SUM, MPI_COMM_WORLD)};

        if (domain.rank() == 0) {
            if (i >= scale_start && i % scale_interval == 0) {
                double strain = total_scaled / whisker_z;
                double stress = (global_force1 - global_force2) / area;
                data << i <<  "," << global_total_energy <<  "," << global_kinetic_energy
                << "," << global_potential_energy << "," << temp <<"," << strain  << "," << stress << std::endl;
            }
            if (i % 10 == 0) write_xyz(traj, atoms);

            // To keep track of which iteration we are in
            if (!(i & (i-1))) std::cout << "i:" << i << std::endl;
        }
        i++;
    }
    data.close();
    traj.close();
    MPI_Finalize();
}
