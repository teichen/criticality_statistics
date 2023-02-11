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
    /* set lattice length under the coarse grain scale, b
    */
    L = L / b;
}

double SimSpace::separation(int* r1, int* r2)
{
    /* calculate a lattice separation between two d=3 positions
    */
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
    /* unpack a flattened position (single integer) into
       three coordinates subject to the lattice length
    */
    r[2] = n % L;
    r[1] = (int)((n - r[2]) / L) % L;
    r[0] = (int)((n - r[2] - r[1] * L) / pow(L, 2)) % L;
}

int SimSpace::flatten_position(int i, int j, int k)
{
    /* flatten a three coordinate position
    */
    int n;
    n = (int)(i * pow(L, 2) + j * L + k);

    return n;
}

void SimSpace::nearest_neighbors(int* r, int nn[3])
{
    /* fill nn with an array of flattened lattice sites for
       the three unique nearest-neighbors
    */
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

void SimSpace::diagonal_in_plane(int* r, int nn[6])
{
    /* fill nn with an array of diagonal in plane neighbors
       organized according to the following:
       sum_{diagonal in plane}\, n_i * n_j
       directions: n_(1,1,0), n_(1,-1,0),
                   n_(1,0,1), n_(1,0,-1),
                   n_(0,1,1), n_(0,1,-1)
       return n_j only
    */
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
    /* fill nn with an array of cubic diagonal neighbors
       organized according to the following:
       sum_{cubic diagonal}\, n_i * n_j
       directions: n_(1,1,1), n_(1,1,-1),
                   n_(1,-1,-1), n_(1,-1,1)
       return n_j only
    */

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
            if (((i>0) && (r[i] == (L-1))) || ((i<0) && (r[i] == 0)))
            {
                nn[k] += (int)(i*(1-L) * L);
            }
            else
            {
                nn[k] += (int)(i * L);
            }
            if (((j>0) && (r[j] == (L-1))) || ((j<0) && (r[j] == 0)))
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
    /* fill nn with an array of principal plane neighbors
       organized according to the following:
       sum_{principal planes}\, n_i * n_j * n_k * n_l
       directions: n_(1,0,0), n_(0,1,0),
                   n_(1,1,0), n_(1,0,0),
                   n_(0,0,1), n_(1,0,1),
                   n_(0,1,0), n_(0,0,1),
                   n_(0,1,1)
    */
    int i, j, k, di, dj;
    for (i=0; i<9; i++)
    {
        nn[i] = (int)(r[0]*pow(L,2) + r[1]*L + r[2]);
    }
    for (i=0; i<2; i++)
    {
        for (j=(i+1); j<3; j++)
        {
            k = 0;
            for (di=0; di<2; di++)
            {
                for (dj=0; dj<2; dj++)
                {
                    if (di>0)
                    {
                        if (r[i] == (L-1))
                        {
                            nn[k] += (int)((1-L) * pow(L, 2-i));
                        }
                        else
                        {
                            nn[k] += (int)(1 * pow(L, 2-i));
                        }
                    }
                    if (dj>0)
                    {
                        if (r[j] == (L-1))
                        {
                            nn[k] += (int)((1-L) * pow(L, 2-j));
                        }
                        else
                        {
                            nn[k] += (int)(1 * pow(L, 2-j));
                        }
                    }
                    k++;
                }
            }
        }
    }
}

void SimSpace::diagonal_planes(int* r, int nn[18])
{
    /* fill nn with an array of diagonal plane neighbors
       organized according to the following:
       sum_{diagonal planes}\, n_i * n_j * n_k * n_l
       directions: n_(1,0,0), n_(0,1,1), n_(1,1,1),
                   n_(1,0,0), n_(0,1,-1), n_(1,1,-1),
                   n_(0,1,0), n_(1,0,1), n_(1,1,1),
                   n_(0,1,0), n_(1,0,-1), n_(1,1,-1),
                   n_(0,0,1), n_(1,1,0), n_(1,1,1),
                   n_(0,0,1), n_(1,-1,0), n_(1,-1,1)
    */
    int i, j, k;
    for (i=0; i<18; i++)
    {
        nn[i] = (int)(r[0]*pow(L,2) + r[1]*L + r[2]);
    }
    // start with the first of each of six triplets
    for (i=0; i<3; i++)
    {
        if (r[i] == (L-1))
        {
            nn[6*i]   += (int)((1-L) * pow(L, 2-i));
            nn[6*i+3] += (int)((1-L) * pow(L, 2-i));
        }
        else
        {
            nn[6*i]   += (int)(1 * pow(L, 2-i));
            nn[6*i+3] += (int)(1 * pow(L, 2-i));
        }
    }
    // now for the second of each triplet
    k = 3*4 + 1;
    for (i=0; i<2; i++)
    {
        for (j=i+1; j<3; j++)
        {
            if (r[i] == (L-1))
            {
                nn[k]   += (int)((1-L) * pow(L, 2-i));
                nn[k+3] += (int)((1-L) * pow(L, 2-i));
            }
            else
            {
                nn[k]   += (int)(1 * pow(L, 2-i));
                nn[k+3] += (int)(1 * pow(L, 2-i));
            }
            if (r[j] == (L-1))
            {
                nn[k]   += (int)((1-L) * pow(L, 2-j));
                nn[k+3] += (int)(-1 * pow(L, 2-j));
            }
            else if (r[j] == 0)
            {
                nn[k]   += (int)(1 * pow(L, 2-j));
                nn[k+3] += (int)((L-1) * pow(L, 2-j));
            }
            else
            {
                nn[k]   += (int)(1 * pow(L, 2-j));
                nn[k+3] += (int)(-1 * pow(L, 2-j));
            }
            k = k - 3*2;
        }
    }
}

void SimSpace::tetrahedral_vertices(int* r, int nn[6])
{
    /* fill nn with tetrahedral vertex neighbors
       organized according to the following:
       sum_{tetrahedral vertices}\, n_i * n_j * n_k * n_l
       directions: n_(1,0,1), n_(1,1,0), n_(0,1,1),
                   n_(1,0,1), n_(1,-1,0), n_(0,-1,1)
    */
    int i, j, k;
    for (i=0; i<6; i++)
    {
        nn[i] = (int)(r[0]*pow(L,2) + r[1]*L + r[2]);
    }
    for (k=0; k<2; k++)
    {
        if (r[0] == (L-1))
        {
            nn[3*k+0] += (int)((1-L) * pow(L, 2));
            nn[3*k+1] += (int)((1-L) * pow(L, 2));
        }
        else
        {
            nn[3*k+0] += (int)(1 * pow(L, 2));
            nn[3*k+1] += (int)(1 * pow(L, 2));
        }
        if (r[2] == (L-1))
        {
            nn[3*k+0] += (int)(1-L);
            nn[3*k+2] += (int)(1-L);
        }
        else
        {
            nn[3*k+0] += 1;
            nn[3*k+2] += 1;
        }
    }
    if (r[1] == (L-1))
    {
        nn[1] += (int)((1-L) * L);
        nn[2] += (int)((1-L) * L);
    }
    else
    {
        nn[1] += 1;
        nn[2] += 1;
    }
    if (r[1] == 0)
    {
        nn[1+3] += (int)((L-1) * L);
        nn[2+3] += (int)((L-1) * L);
    }
    else
    {
        nn[1+3] += -1;
        nn[2+3] += -1;
    }
}

void SimSpace::next_nearest_neighbors(int* r, int nn[3])
{
    /* fill nn with an array of next nearest neighbors
       organized according to the following:
       sum_NNN\, n_i * n_i+2 (i=x,y,z)
    */
    int i;
    for (i=0; i<dim; i++)
    {
        nn[i] = (int)(r[0]*pow(L,2) + r[1]*L + r[2]);
        if (r[i] >= (L-2))
        {
            nn[i] += (int)((2-L) * pow(L, 2-i));
        }
        else
        {
            nn[i] += (int)(2 * pow(L, 2-i));
        }
    }
}

SimSpace::~SimSpace()
{
}
