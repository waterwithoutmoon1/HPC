//
// Created by 13396 on 24-9-1.
//

#ifndef IH_H
#define IH_H

#include <vector>
#include "vector.h"

struct Particle
{
    char s[10];
    Vector v;
}*p;

struct Edge
{
    int i;
    int j;
}e[30];

struct Facet
{
    int i;
    int j;
    int k;
}f[30];

std::vector<Particle> generate_cluster(int n);

#endif
