//
// Created by 13396 on 2024/6/13.
//

#include <iomanip>
#include <iostream>

#include <filesystem>

#include "xyz.h"
#include <sstream>

std::tuple<Names_t, Positions_t> read_xyz(const std::string &filename) {
    std::ifstream file(filename);

    if (file.is_open()) {
        std::string line;

        // First line contains number of atoms
        std::getline(file, line);
        int nb_atoms;
        std::istringstream(line) >> nb_atoms;

        // Second line contains a comment - in extended XYZ this line contains auxiliary information
        std::getline(file, line);

        // Data structures for names and positions
        Names_t names(nb_atoms);
        Eigen::Array3Xd positions(3, nb_atoms);
        positions.setZero();

        // Now follows a line for each atom
        for (int i = 0; i < nb_atoms; ++i) {
            std::getline(file, line);
            std::istringstream(line) >> names[i] >> positions(0, i) >> positions(1, i) >> positions(2, i);
        }

        // Close file, we're done
        file.close();

        return {names, positions};
    } else {
        throw std::runtime_error("Could not open file");
    }
}


std::tuple<Names_t, Positions_t, Velocities_t> read_xyz_with_velocities(const std::string &filename) {
    std::ifstream file(filename);

    if (file.is_open()) {
        std::string line;

        // First line contains number of atoms
        std::getline(file, line);
        int nb_atoms;
        std::istringstream(line) >> nb_atoms;

        // Second line contains a comment - in extended XYZ this line contains
        // auxiliary information
        std::getline(file, line);

        // Data structures for names and positions
        Names_t names(nb_atoms);
        Eigen::Array3Xd positions(3, nb_atoms);
        Eigen::Array3Xd velocities(3, nb_atoms);
        positions.setZero();
        velocities.setZero();

        // Now follows a line for each atom
        for (int i = 0; i < nb_atoms; ++i) {
            std::getline(file, line);
            std::istringstream(line) >> names[i] >> positions(0, i) >> positions(1, i) >> positions(2, i)
                                     >> velocities(0, i) >> velocities(1, i) >> velocities(2, i);
        }

        // Close file, we're done
        file.close();

        return {names, positions, velocities};
    } else {
        std::cerr << "Could not open file: " << filename << std::endl;

        // Print the current working directory
        std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;

        throw std::runtime_error("Could not open file");
    }
}


void write_xyz(std::ofstream &file, Atoms &atoms) {
    // Number of atoms
    file << atoms.nb_atoms() << std::endl;

    // Comment line
    file << std::endl;

    // Element name, position
    for (int i = 0; i < atoms.nb_atoms(); ++i) {
        file << std::setw(2) << "Au" << " "
             << std::setw(10) << atoms.positions.col(i).transpose()
             << std::setw(10) << atoms.velocities.col(i).transpose()
             << std::endl;
    }
}


void write_xyz(const std::string &filename, Atoms &atoms) {
    // Open file for writing
    std::ofstream file(filename);

    // Write XYZ
    write_xyz(file, atoms);

    // Close file, we're done
    file.close();
}