//
// Created by 13396 on 24-9-1.
//
#include "ih.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

// Function to initialize basic vectors and edges
void init_vectors_and_edges(Vector b[], Edge e[], Facet f[]) {
  p = NULL;

  //
  // Initialize basic vectors.
  //
  const double HT = ( sqrt(5.0) + 1.0 ) / 4.0;   // half Tao

  b[0] = Vector( HT, 0.0, 0.5 );
  b[1] = Vector( HT, 0.0, -0.5 );
  b[2] = Vector( 0.5, HT, 0.0 );
  b[3] = Vector( -0.5, HT, 0.0 );
  b[4] = Vector( 0.0, 0.5, HT );
  b[5] = Vector( 0.0, -0.5, HT );
  b[6] = Vector( 0.5, -HT, 0.0 );
  b[7] = Vector( 0.0, 0.5, -HT );
  b[8] = Vector( -HT, 0.0, 0.5 );
  b[9] = Vector( 0.0, -0.5, -HT );
  b[10] = Vector( -HT, 0.0, -0.5 );
  b[11] = Vector( -0.5, -HT, 0.0 );

  //
  // Initialize 30 edges
  //
  e[0].i = 0;  e[0].j = 1;
  e[1].i = 0;  e[1].j = 2;
  e[2].i = 0;  e[2].j = 4;
  e[3].i = 0;  e[3].j = 5;
  e[4].i = 0;  e[4].j = 6;

  e[5].i = 10;  e[5].j = 3;
  e[6].i = 10;  e[6].j = 7;
  e[7].i = 10;  e[7].j = 8;
  e[8].i = 10;  e[8].j = 9;
  e[9].i = 10;  e[9].j = 11;

  e[10].i = 1;  e[10].j = 2;
  e[11].i = 1;  e[11].j = 6;
  e[12].i = 1;  e[12].j = 7;
  e[13].i = 1;  e[13].j = 9;

  e[14].i = 8;  e[14].j = 3;
  e[15].i = 8;  e[15].j = 4;
  e[16].i = 8;  e[16].j = 5;
  e[17].i = 8;  e[17].j = 11;

  e[18].i = 2;  e[18].j = 3;
  e[19].i = 2;  e[19].j = 4;
  e[20].i = 2;  e[20].j = 7;

  e[21].i = 11;  e[21].j = 5;
  e[22].i = 11;  e[22].j = 6;
  e[23].i = 11;  e[23].j = 9;

  e[24].i = 6;  e[24].j = 5;
  e[25].i = 6;  e[25].j = 9;

  e[26].i = 3;  e[26].j = 4;
  e[27].i = 3;  e[27].j = 7;

  e[28].i = 7;  e[28].j = 9;

  e[29].i = 5;  e[29].j = 4;

  //
  // Initialize 20 facets
  //
  f[0].i = 0;   f[0].j = 1;   f[0].k = 2;
  f[1].i = 0;   f[1].j = 2;   f[1].k = 4;
  f[2].i = 0;   f[2].j = 4;   f[2].k = 5;
  f[3].i = 0;   f[3].j = 5;   f[3].k = 6;
  f[4].i = 0;   f[4].j = 1;   f[4].k = 6;

  f[5].i = 10;   f[5].j = 3;   f[5].k = 7;
  f[6].i = 10;   f[6].j = 3;   f[6].k = 8;
  f[7].i = 10;   f[7].j = 8;   f[7].k = 11;
  f[8].i = 10;   f[8].j = 9;   f[8].k = 11;
  f[9].i = 10;   f[9].j = 7;   f[9].k = 9;

  f[10].i = 1;   f[10].j = 2;   f[10].k = 7;
  f[11].i = 1;   f[11].j = 7;   f[11].k = 9;
  f[12].i = 1;   f[12].j = 6;   f[12].k = 9;

  f[13].i = 8;   f[13].j = 5;   f[13].k = 11;
  f[14].i = 8;   f[14].j = 4;   f[14].k = 5;
  f[15].i = 8;   f[15].j = 3;   f[15].k = 4;

  f[16].i = 2;   f[16].j = 3;   f[16].k = 7;
  f[17].i = 2;   f[17].j = 3;   f[17].k = 4;

  f[18].i = 11;   f[18].j = 5;   f[18].k = 6;
  f[19].i = 11;   f[19].j = 6;   f[19].k = 9;
}

// Function to calculate the number of particles on the nth layer
int np( int n )
{
    if( n<0 ) return -1;
    else if( n==0 ) return 1;
    else if( n==1 ) return 12;
    else if( n==2 ) return 42;
    else
    {
        int count = 0;

        count += 12;   // edge particles
        count += (n-1)*30;   // side particles
        for( int i=1; i<=n-2; i++ ) count += i*20;   // body particles
        return count;
    }
}

std::vector<Particle> generate_cluster(int n) {
    // Generate particles as per the logic in `ih.C`
    std::vector<Particle> particles;

    if (n < 0) return particles;

    int count = 0;

    // Initialize the basic vectors and edges
    Vector b[12];
    Edge e[30];
    Facet f[20];
    init_vectors_and_edges(b, e, f);

    int num_particles = np(n);
    particles.reserve(num_particles);

    if (n == 0) {
        Particle center;
        strcpy(center.s, "H");
        center.v = Vector(0.0, 0.0, 0.0);
        particles.push_back(center);
        return particles;
    }

    // Generate edge particles
    for (int i = 0; i < 12; ++i) {
        Particle particle;
        strcpy(particle.s, "Ag");
        particle.v = b[i] * n;
        particles.push_back(particle);
    }

    // Generate side particles
    if (n >= 2) {
        for (int i = 0; i < 30; ++i) {
            Vector e1 = b[e[i].i] * n;
            Vector e2 = b[e[i].j] * n;
            for (int j = 1; j <= n - 1; ++j) {
                Particle particle;
                strcpy(particle.s, "Cd");
                particle.v = e1 + (e2 - e1) * j / double(n);
                particles.push_back(particle);
            }
        }
    }

    // Generate body particles
    if (n >= 3) {
        for (int i = 0; i < 20; ++i) {
            Vector e1 = b[f[i].i] * n;
            Vector e2 = b[f[i].j] * n;
            Vector e3 = b[f[i].k] * n;
            for (int j = 1; j <= n - 2; ++j) {
                Vector v1 = e1 + (e2 - e1) * (j + 1) / double(n);
                Vector v2 = e1 + (e3 - e1) * (j + 1) / double(n);
                for (int k = 1; k <= j; ++k) {
                    Particle particle;
                    strcpy(particle.s, "Au");
                    particle.v = v1 + (v2 - v1) * k / double(j + 1);
                    particles.push_back(particle);
                }
            }
        }
    }

    return particles;
}