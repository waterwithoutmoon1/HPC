//
// Created by 13396 on 2024/6/19.
//

#ifndef LJ_NEIGHBORLIST_H
#define LJ_NEIGHBORLIST_H
#include <atoms.h>
#include <neighbors.h>

double lj_neighborlist(Atoms &atoms, NeighborList &neighbor_list, double epsilon = 1.0, double sigma = 1.0, double cutoff = 0);

#endif //LJ_NEIGHBORLIST_H
