//
// Created by 13396 on 2024/6/5.
//

#ifndef ATOMS_H
#define ATOMS_H
#include <Eigen/Dense>

using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
using Masses_t = Eigen::VectorXd;;
using Names_t = std::vector<std::string>;

class Atoms {
public:
    Eigen::Array<bool, Eigen::Dynamic, 1> is_ghost;
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Masses_t masses;
    Names_t names;
    Eigen::ArrayXd per_atom_potential_energy;

    Atoms(const Positions_t &p) :
            positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{p.cols()} {
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
        per_atom_potential_energy.setZero();
        initialize_ghost_flags();
    }

    Atoms(const Positions_t &p, const Velocities_t &v) :
            positions{p}, velocities{v}, forces{3, p.cols()}, masses{p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
        masses.setOnes();
        per_atom_potential_energy.setZero();
        initialize_ghost_flags();
    }

    Atoms(const Names_t &n, const Positions_t &p) :
            names{n},positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{p.cols()} {
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
        per_atom_potential_energy.setZero();
        initialize_ghost_flags();
    }

    [[nodiscard]] size_t nb_atoms() const {
        return positions.cols();
    }

    void resize(const size_t new_size) {
        size_t old_size = nb_atoms();
        positions.conservativeResize(Eigen::NoChange, new_size);
        velocities.conservativeResize(Eigen::NoChange, new_size);
        forces.conservativeResize(Eigen::NoChange, new_size);
        masses.conservativeResize(new_size);
        // Initialize newly added atoms as non-ghosts
        is_ghost.conservativeResize(new_size);
        per_atom_potential_energy.conservativeResize(new_size);
        if (new_size > old_size) {
            is_ghost.segment(old_size, new_size - old_size).setConstant(false);
        }
    }

    // Method to initialize all atoms as non-ghost
    void initialize_ghost_flags() {
        is_ghost.setConstant(nb_atoms(), false);
    }

    bool is_ghost_atom(size_t index) const {
        return is_ghost(index);
    }
};

#endif //ATOMS_H
