//
// Created by 13396 on 2024/5/25.
//
#include <iostream>
#include <lj_direction_summation.cpp>
#include <string>
#include <xyz.h>
#include <chrono>
#include <verlet_2.h>
#include <berendsen_thermostat.h>


void velocity_verlet(Atoms &atoms, double timestep, double epsilon, double sigma, double mass) {
    size_t num_atoms = atoms.nb_atoms();
    // Compute initial forces and potential energy
    double potential_energy = lj_direct_summation(atoms, epsilon, sigma);

    // Update positions
    verlet2_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);

    // Compute forces and potential energy
    potential_energy = lj_direct_summation(atoms, epsilon, sigma);

    // Update velocities
    verlet2_step2(atoms.velocities, atoms.forces, timestep, mass);
}


// Function to run the simulation
double run_simulation(Atoms &atoms, double kB, double epsilon, double sigma, double mass,
                    double target_temp, double tau_T, double timestep, double total_time,
                    const std::string &output_file) {
    int steps = static_cast<int>(total_time / timestep);

    // File to record data
    std::ofstream data(output_file);
    data << "step,total_energy,kinetic_energy,potential_energy,temperature\n";

    auto start = std::chrono::high_resolution_clock::now();

    // Main MD loop
    for (int step = 0; step <= steps; ++step) {
        // Compute potential energy and forces
        double potential_energy = lj_direct_summation(atoms, epsilon, sigma);

        // Compute kinetic energy and temperature
        double kinetic_energy = compute_kinetic_energy(atoms.velocities,mass);
        double total_energy = potential_energy + kinetic_energy;
        double temperature = (2.0 * kinetic_energy) / (3.0 * atoms.nb_atoms() * kB);

        // Record data
        data << step << "," << total_energy << "," << kinetic_energy << ","
             << potential_energy << "," << temperature << "\n";

        // Velocity Verlet integration
        velocity_verlet(atoms, timestep, epsilon, sigma, mass);

        // Apply Berendsen thermostat
        apply_berendsen_thermostat(atoms, kB, mass, target_temp, tau_T, timestep);
    }

    //Caculate the simulation time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Number of atoms: " << atoms.nb_atoms() << ", Simulation time: " << duration.count() << " seconds" << std::endl;
    double elapsed_time = duration.count();

   //Close the write of data
    data.close();

    return elapsed_time;
}

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

int main() {
try{
    //Parameters of atoms
    Positions_t positions;
    Velocities_t velocities;
    Atoms atoms(positions, velocities);
    double epsilon = 1.0;
    double sigma = 1.0;
    double mass = 1.0;
    double kB = 1.0;

    //Parameters of time
    double timestep = 0.001;
    double total_time = 2.5;
    double target_temp = 1;  // Target temperature for the thermostat
    double tau_T = 1.0;  // Thermostat coupling constant

    //Marix of amount of atoms and simulation time
    std::vector<int> atom_counts = {8, 27, 64, 125, 216, 343, 512, 729, 1000};
    std::vector<double> times;

    //Simulation with different numbers of atoms
    for (int num_atoms_per_side = 2; num_atoms_per_side <= 10; ++num_atoms_per_side) {
        double lattice_constant = 1.0;
        create_cubic_lattice(atoms, num_atoms_per_side, lattice_constant);

        std::string output_file = "simulation_data_" + std::to_string(atoms.nb_atoms()) + ".csv";

        double elapsed_time = run_simulation(atoms, kB, epsilon, sigma, mass,
                            target_temp, tau_T, timestep,  total_time,
                            output_file);

        times.push_back(elapsed_time);
    }

    // Save the results to a file
    std::ofstream result_file("simulation_times.csv");
    result_file << "num_atoms,elapsed_time\n";
    for (int i = 0; i < atom_counts.size(); ++i) {
        result_file << atom_counts[i] << "," << times[i] << "\n";
    }
    result_file.close();
} catch (const std::runtime_error &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}