/*
    The blocking class calculates blocking statistics of a 
    binary discretized field.
*/
#include "blocking.h"
#include <cstdlib>
#include <stdlib.h>
#include <fstream>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <sstream>
using std::ifstream;
#include <cstring>

#include <sys/time.h>
#include <unistd.h>

#include <time.h>
#include <math.h>
#include <cmath>

#include <omp.h>

const double PI = 3.1415926535897932384626433832795028841971693;

using namespace std;

blocking::blocking()
{
    lattice.set_dimensions(1);
    L   = lattice.L; // length of lattice (number of sites)
    dim = lattice.dim; // dimensionality of lattice
    nL  = pow(L,dim);
}

void blocking::coarse_grain_field(int* n, int b_blocking)
{
    // TODO: call from Executive which loads configs from ./netlib/n*.dat
    b = b_blocking;

    coarse_lattice.set_dimensions(b);
    Lb = coarse_lattice.L;

    mem_test = false;
    initarrays();

    int ub[nL];
    lattice_map(bshift, ub); // the coarse cell that each original lattice site belongs to

    std::ostringstream filenameStream;
    std:string filename;

    ifstream fin;
    string line;

    std::ofstream ndump;

    int i,j,k;
    int ib,jb,kb;

    // collective variables for each MC config
    // n_collective[0] = sum\, n_i
    // n_collective[1] = sum_NN\, n_i * n_i+1 (i=x,y,z)
    // n_collective[2] = sum_{diagonal in plane}\, n_i * n_j
    // n_collective[3] = sum_{cubic diagonal}\, n_i * n_j
    // n_collective[4] = sum_{principal planes}\, n_i * n_j * n_k * n_l
    // n_collective[5] = sum_{diagonal planes}\, n_i * n_j * n_k * n_l
    // n_collective[6] = sum_{tetrahedral vertices}\, n_i * n_j * n_k * n_l
    // n_collective[7] = sum_NNN\, n_i * n_i+2 (i=x,y,z)

    int celli,cellj,cellk,celll;

    double n_tmp;

    int nn[18];
    int nn_vals[18];
    int ri[3];

    double n_i,n_j,n_k,n_l;

    if (b>1)
    {
        bshift = rand() % b;
    }
    else
    {
        bshift = 0;
    }

    i = 0; j = 0; k = 0;

    ib = 0; jb = 0; kb = 0;

    for (i=0; i<(int)(pow(Lb,dim)); i++)
    {
        nb[i] = 0.0;
    }

    double n_1, n0, n1, n2, n3, n4, n5, n6;
    n_1 = 0.0; n0 = 0.0; n1 = 0.0; n2 = 0.0; n3 = 0.0; n4 = 0.0; n5 = 0.0; n6 = 0.0;

    n_tmp = 0.0;

    n_i = 0.0; n_j = 0.0; n_k = 0.0; n_l = 0.0;

    #pragma omp parallel for shared (n,nb,ub) private (i,n_tmp)
    for (i=0; i<(int)(pow(L,dim)); i++)
    {
    #pragma omp atomic read
        n_tmp = n[i];
    #pragma omp atomic write
        nb[ub[i]] += n_tmp;
    }

    #pragma omp parallel for shared (nb) private (i,n_tmp) reduction(+:n_1)
    for (i=0; i<(int)(pow(Lb,dim)); i++)
    {
    #pragma omp atomic read
        n_tmp = nb[i];
    #pragma omp atomic write
        nb[i] = n_tmp / (double)(pow(b,dim));

        n_1 += n_tmp / (double)(pow(b,dim));
    }

    // TODO: write out coarse grained field to nb.dat

    #pragma omp parallel for shared (nb) private (i,j,ri,nn,nn_vals,n_i,n_j) reduction(+:n0)
    for (i=0; i<(int)(pow(Lb,dim)); i++)
    {
        // sum_NN\, n_i * n_i+1 (i=x,y,z)
        lattice.unpack_position(i, ri);
        lattice.nearest_neighbors(ri, nn);
    
        for (j=0; j<3; j++)
        {
    #pragma omp atomic read
            nn_vals[j] = nb[nn[j]];
        }

    #pragma omp atomic read
        n_i = nb[i];
        for (j=0; j<3; j++)
        {
            n_j = nn_vals[j];

            n0 += n_i * n_j;
        }
    }

    #pragma omp parallel for shared (nb) private (i,j,ri,nn,nn_vals,n_i,n_j) reduction(+:n1)
    for (i=0; i<(int)(pow(Lb,dim)); i++)
    {
        // sum_{diagonal in plane}\, n_i * n_j

        lattice.unpack_position(i, ri);
        lattice.diagonal_in_plane(ri, nn);
    
        for (j=0; j<6; j++)
        {
    #pragma omp atomic read
            nn_vals[j] = nb[nn[j]];
        }

    #pragma omp atomic read
        n_i = nb[i];
        for (j=0; j<6; j++)
        {
            n_j = nn_vals[j];

            n1 += n_i * n_j;
        }
    }

    #pragma omp parallel for shared (nb) private (i,j,ri,nn,nn_vals,n_i,n_j) reduction(+:n2)
    for (i=0; i<(int)(pow(Lb,dim)); i++)
    {
        // sum_{cubic diagonal}\, n_i * n_j
        lattice.unpack_position(i, ri);
        lattice.cubic_diagonal(ri, nn);
    
        for (j=0; j<4; j++)
        {
    #pragma omp atomic read
            nn_vals[j] = nb[nn[j]];
        }

    #pragma omp atomic read
        n_i = nb[i];
        for (j=0; j<4; j++)
        {
            n_j = nn_vals[j];

            n2 += n_i * n_j;
        }
    }

    #pragma omp parallel for shared (nb) private (i,j,ri,nn,nn_vals,n_i,n_j,n_k,n_l) reduction(+:n3)
    for (i=0; i<(int)(pow(Lb,dim)); i++)
    {
        // sum_{principal planes}\, n_i * n_j * n_k * n_l

        lattice.unpack_position(i, ri);
        lattice.principal_planes(ri, nn);
    
        for (j=0; j<9; j++)
        {
    #pragma omp atomic read
            nn_vals[j] = nb[nn[j]];
        }

    #pragma omp atomic read
        n_i = nb[i];
        for (j=0; j<3; j++) // three triplets
        {
            n_j = nn_vals[3*j+0];
            n_k = nn_vals[3*j+1];
            n_l = nn_vals[3*j+2];

            n3 += n_i * n_j * n_k * n_l;
        }
    }

    #pragma omp parallel for shared (nb) private (i,j,ri,nn,nn_vals,n_i,n_j,n_k,n_l) reduction(+:n4)
    for (i=0; i<(int)(pow(Lb,dim)); i++)
    {
        // sum_{diagonal planes}\, n_i * n_j * n_k * n_l
        lattice.unpack_position(i, ri);
        lattice.diagonal_planes(ri, nn);
    
        for (j=0; j<18; j++)
        {
    #pragma omp atomic read
            nn_vals[j] = nb[nn[j]];
        }

    #pragma omp atomic read
        n_i = nb[i];
        for (j=0; j<6; j++) // three triplets
        {
            n_j = nn_vals[3*j+0];
            n_k = nn_vals[3*j+1];
            n_l = nn_vals[3*j+2];

            n4 += n_i * n_j * n_k * n_l;
        }
    }

    #pragma omp parallel for shared (nb) private (i,j,ri,nn,nn_vals,n_i,n_j,n_k,n_l) reduction(+:n5)
    for (i=0; i<(int)(pow(Lb,dim)); i++)
    {
        // sum_{tetrahedral vertices}\, n_i * n_j * n_k * n_l
        lattice.unpack_position(i, ri);
        lattice.tetrahedral_vertices(ri, nn);
    
        for (j=0; j<6; j++)
        {
    #pragma omp atomic read
            nn_vals[j] = nb[nn[j]];
        }

    #pragma omp atomic read
        n_i = nb[i];
        for (j=0; j<2; j++) // three triplets
        {
            n_j = nn_vals[3*j+0];
            n_k = nn_vals[3*j+1];
            n_l = nn_vals[3*j+2];

            n5 += n_i * n_j * n_k * n_l;
        }
    }

    #pragma omp parallel for shared (nb) private (i,j,ri,nn,nn_vals,n_i,n_j) reduction(+:n6)
    for (i=0; i<(int)(pow(Lb,dim)); i++)
    {
        lattice.unpack_position(i, ri);
        lattice.next_nearest_neighbors(ri, nn);
    
        for (j=0; j<3; j++)
        {
    #pragma omp atomic read
            nn_vals[j] = nb[nn[j]];
        }

    #pragma omp atomic read
        n_i = nb[i];
        for (j=0; j<3; j++)
        {
            n_j = nn_vals[j];

            n6 += n_i * n_j;
        }
    }

    n_collective[0] = n_1;
    n_collective[1] = n0;
    n_collective[2] = n1;
    n_collective[3] = n2;
    n_collective[4] = n3;
    n_collective[5] = n4;
    n_collective[6] = n5;
    n_collective[7] = n6;

    // TODO: save results to disk, write n_collective to nvec.dat

}

void blocking::lattice_map(int bshift, int* ub)
{
    int ri[3];

    int i, j;
    for (i=0; i<(int)(pow(L,dim)); i++)
    {
        lattice.unpack_position(i, ri);

        // randomized sampling, shift the coarse grid boundaries (isotropic)
        ri[0] += bshift;
        ri[1] += bshift;
        ri[2] += bshift;

        // periodic boundary conditions
        for (j=0; j<dim; j++)
        {
            if (ri[j] >= L)
            {
                ri[j] -= L;
            }
        }

        ri[2] /= b;
        ri[1] /= b;
        ri[0] /= b;

        ub[i] = (int)(ri[0] * pow(Lb, 2) + ri[1] * Lb + ri[2]);
    }
}

void blocking::initarrays()
{
    nb = (double*) calloc (pow(Lb, dim), sizeof(double));

    mem_test = true;
}

blocking::~blocking()
{
    if(mem_test==true)
    {
    delete [] nb;
    }
}

