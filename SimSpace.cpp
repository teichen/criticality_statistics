/*
    The SimSpace class sets the dimensions of the
    simulation volume.
*/
#include "SimSpace.h"
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <math.h>
#include <iostream>

using namespace std;

SimSpace::SimSpace()
{
    L   = 8;
    dim = 3;
}

void SimSpace::set_dimensions(int b)
{
    L = L / b;
}

double SimSpace::separation(int* r1, int* r2)
{
    double d[3];
    d[0] = 0.0; d[1] = 0.0; d[2] = 0.0;

    double r;
    int dimR;

    for (dimR = 0; dimR < dim; dimR++)
    {
        if ((double)(r1[dimR] - r2[dimR]) > (double)(L/2))
        {
            d[dimR] = (double)(r1[dimR] - r2[dimR] - L);
        }
        else if((double)(r1[dimR] - r2[dimR]) < (double)(-L/2))
        {
            d[dimR] = (double)(r1[dimR] - r2[dimR] + L);
        }
        else
        {
            d[dimR] = (double)(r1[dimR] - r2[dimR]);
        }
    }

    r = sqrt(pow(d[0], 2) + pow(d[1], 2) + pow(d[2], 2));

    return r;
}

void SimSpace::unpack_position(int n, int r[3])
{
    r[2] = n % L;
    r[1] = (int)((n - r[2]) / L) % L;
    r[0] = (int)((n - r[2] - r[1] * L) / pow(L, 2)) % L;
}

int SimSpace::flatten_position(int i, int j, int k)
{
    // position n -> (x*L^2+y*L+z)
    int n;
    n = (int)(i * pow(L, 2) + j * L + k);

    return n;
}

void SimSpace::nearest_neighbors(int* r, int nn[3])
{
    int i;
    for (i=0; i<dim; i++)
    {
        nn[i] = (int)(r[0]*pow(L,2) + r[1]*L + r[2]);
        if (r[i] == (L-1))
        {
            nn[i] += (int)((1-L) * pow(L, 2-i));
        }
        else
        {
            nn[i] += (int)(1 * pow(L, 2-i));
        }
    }
}

void SimSpace::nearest_neighbor_values(int* n, int* nn, int nn_vals[3])
{
    nn_vals[0] = n[nn[0]];
    nn_vals[1] = n[nn[1]];
    nn_vals[2] = n[nn[2]];
}

void SimSpace::diagonal_in_plane(int* r, int nn[6])
{
    // sum_{diagonal in plane}\, n_i * n_j
    // directions: n_(1,1,0), n_(1,-1,0),
    //             n_(1,0,1), n_(1,0,-1),
    //             n_(0,1,1), n_(0,1,-1)
    // return n_j only

    int i, j;
    for (i=0; i<6; i++)
    {
        nn[i] = (int)(r[0]*pow(L,2) + r[1]*L + r[2]);
    }
    for (i=0; i<2; i++)
    {
        for (j=(i+1); j<3; j++)
        {
            if (r[i] == (L-1))
            {
                nn[i] += (int)((1-L) * pow(L, 2-i));
            }
            else
            {
                nn[i] += (int)(1 * pow(L, 2-i));
            }
            if (r[j] == 0)
            {
                nn[j] += (int)((L-1) * pow(L, 2-j));
            }
            else
            {
                nn[j] += (int)(-1 * pow(L, 2-j));
            }
        }
    }
}

void SimSpace::cubic_diagonal(int* r, int nn[4])
{
    // sum_{cubic diagonal}\, n_i * n_j
    // directions: n_(1,1,1), n_(1,1,-1),
    //             n_(1,-1,-1), n_(1,-1,1)
    // return n_j only

    int i, j, k;
    for (i=0; i<4; i++)
    {
        nn[i] = (int)(r[0]*pow(L,2) + r[1]*L + r[2]);
        if (r[0] == (L-1))
        {
            nn[i] += (int)((1-L) * pow(L, 2));
        }
        else
        {
            nn[i] += (int)(1 * pow(L, 2));
        }
    }

    k = 0;
    for (i=-1; i<=1; i=i+2)
    {
        for (j=-1; j<=1; j=j+2)
        {
            if ((i>0) && (r[i] == (L-1))) || ((i<0) && (r[i] == 0))
            {
                nn[k] += (int)(i*(1-L) * L);
            }
            else
            {
                nn[k] += (int)(i * L);
            }
            if ((j>0) && (r[j] == (L-1))) || ((j<0) && (r[j] == 0))
            {
                nn[k] += (int)(j*(1-L));
            }
            else
            {
                nn[k] += (int)(j);
            }
            k++;
        }
    }
}

void SimSpace::principal_planes(int* r, int nn[9])
{
    // sum_{principal planes}\, n_i * n_j * n_k * n_l
    // e.g. n_(0,0,0) * n_(0,1,0) * n_(1,0,0) * n_(1,1,0)
n_(1,0,0)
n_(0,1,0)
n_(1,1,0)
n_(1,0,0)
n_(0,0,1)
n_(1,0,1)
n_(0,1,0)
n_(0,0,1)
n_(0,1,1)
}

void SimSpace::diagonal_planes(int* r, int nn[18])
{
    // sum_{diagonal planes}\, n_i * n_j * n_k * n_l
    // e.g. n_(0,0,0) * n_(1,0,0) * n_(0,1,1) * n_(1,1,1)
n_(1,0,0)
n_(0,1,1)
n_(1,1,1)
n_(1,0,0)
n_(0,1,-1)
n_(1,1,-1)
n_(0,1,0)
n_(1,0,1)
n_(1,1,1)
n_(0,1,0)
n_(1,0,-1)
n_(1,1,-1)
n_(0,0,1)
n_(1,1,0)
n_(1,1,1)
n_(0,0,1)
n_(1,-1,0)
n_(1,-1,1)
}

void SimSpace::tetrahedral_vertices(int* r, int nn[6])
{
    // sum_{tetrahedral vertices}\, n_i * n_j * n_k * n_l
    // e.g. n_(0,0,0) * n_(1,0,1) * n_(0,1,1) * n_(1,1,0)
n_(1,0,1)
n_(0,1,1)
n_(1,1,0)
n_(1,0,1)
n_(0,-1,1)
n_(1,-1,0)
}

void SimSpace::next_nearest_neighbors(int* r, int nn[6])
{
    // sum_NNN\, n_i * n_i+2 (i=x,y,z)
}

SimSpace::~SimSpace()
{
}
