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


blocking::blocking(int b_blocking)
{
    b = b_blocking;

    coarse_lattice.set_dimensions(b);
    Lb = coarse_lattice.L;

    int ub[nL];
    lattice_map(bshift, ub); // the coarse cell that each original lattice site belongs to

    int n_configs;
    n_configs = 20;

    std::ostringstream filenameStream;
    std:string filename;

    ifstream fin;
    string line;

    std::ofstream ndump;

    int i,j,k;
    int ib,jb,kb;

    // read in data

    int n[nL];

    // coarse-grained n (average per cell)
    double nb[(int)pow(Lb,dim)];

    // collective variables for each MC config
    double n_1; // sum\, n_i
    double n0;  // sum_NN\, n_i * n_i+1 (i=x,y,z)
    double n1;  // sum_{diagonal in plane}\, n_i * n_j
    // e.g. n_(0,0,0) * n_(1,1,0)
    double n2; // sum_{cubic diagonal}\, n_i * n_j
    // e.g. n_(0,0,0) * n_(1,1,1)
    double n3; // sum_{principal planes}\, n_i * n_j * n_k * n_l
    // e.g. n_(0,0,0) * n_(0,1,0) * n_(1,0,0) * n_(1,1,0)
    double n4; // sum_{diagonal planes}\, n_i * n_j * n_k * n_l
    // e.g. n_(0,0,0) * n_(1,0,0) * n_(0,1,1) * n_(1,1,1)
    double n5; // sum_{tetrahedral vertices}\, n_i * n_j * n_k * n_l
    // e.g. n_(0,0,0) * n_(1,0,1) * n_(0,1,1) * n_(1,1,0)
    double n6; // sum_NNN\, n_i * n_i+2 (i=x,y,z)

    double n_1vec[n_configs];
    double n0vec[n_configs];
    double n1vec[n_configs];
    double n2vec[n_configs];
    double n3vec[n_configs];
    double n4vec[n_configs];
    double n5vec[n_configs];
    double n6vec[n_configs];

    int celli,cellj,cellk,celll;

    double n_tmp;

    int nn[18];
    int nn_vals[18];
    int ri[3];

    double nn_i,nn_j,nn_k,nn_l;

    int config;

    for (config=0; config<n_configs; config++)
    {
        if (b>1)
        {
            bshift = rand() % b;
        }
        else
        {
            bshift = 0;
        }

        // read in data

        filenameStream << "./netlib/n" << config << ".dat";
        filename = filenameStream.str();

        fin.open(filename.c_str());

        j = 0;
        while( getline(fin,line) )
        {
            n[j] = atoi (line.c_str());
            j++;
        }

        fin.close();
        filenameStream.str("");

        i = 0; j = 0; k = 0;

        ib = 0; jb = 0; kb = 0;

        for (i=0; i<(int)(pow(Lb,dim)); i++)
        {
            nb[i] = 0.0;
        }

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

            n_1 += nn/(double)(pow(b,dim));
        }

        ndump.open("nb.dat", std::ios_base::app);

        for (i=0; i<(int)(pow(Lb,dim)); i++)
        {
            ndump << nb[i] << "\n";
        }
        ndump.close();

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

    #pragma omp parallel for shared (nb) private (i,j,ri,nn,nn_vals,nn_i,nn_j,nn_k,nn_l) reduction(+:n3)
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

    #pragma omp parallel for shared (nb) private (i,j,ri,nn,nn_vals,nn_i,nn_j,nn_k,nn_l) reduction(+:n4)
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

    #pragma omp parallel for shared (nb) private (i,j,ri,nn,nn_vals,nn_i,nn_j,nn_k,nn_l) reduction(+:n5)
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

    #pragma omp parallel for shared (nb) private (i,j,k,celli,cellj,nn_i,nn_j) reduction(+:n6)
        for (i=0; i<Lb; i++)
        {
    void next_nearest_neighbors(int*, int[6]);


            for (j=0; j<Lb; j++)
            {
                for (k=0; k<Lb; k++)
                {
                    // sum_NNN\, n_i * n_i+2 (i=x,y,z)

                    celli = (int)(i*pow(Lb,2)+j*(Lb)+k);
        #pragma omp atomic read
                    nn_i = nb[celli];

                    //cellj = get_cell(i,j,k,2,0,0,b); 
                    if (i==((Lb)-1))
                    {
                        cellj = get_cell(0,j,k,1,0,0,b);
                    }
                    else
                    {
                        cellj = get_cell(i+1,j,k,1,0,0,b);
                    }

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n6 += nn_i*nn_j; 

                    //cellj = get_cell(i,j,k,0,2,0,b); 
                    if (j==((Lb)-1))
                    {
                        cellj = get_cell(i,0,k,0,1,0,b);
                    }
                    else
                    {
                        cellj = get_cell(i,j+1,k,0,1,0,b);
                    }

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n6 += nn_i*nn_j; 

                    //cellj = get_cell(i,j,k,0,0,2,b); 
                    if (k==((Lb)-1))
                    {
                        cellj = get_cell(i,j,0,0,0,1,b);
                    }
                    else
                    {
                        cellj = get_cell(i,j,k+1,0,0,1,b);
                    }

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n6 += nn_i*nn_j; 
                }
            }
        }

        n_1vec[config] = n_1;
        n0vec[config] = n0;
        n1vec[config] = n1;
        n2vec[config] = n2;
        n3vec[config] = n3;
        n4vec[config] = n4;
        n5vec[config] = n5;
        n6vec[config] = n6;

    }

    // save results to disk

    ndump.open("nvec.dat", std::ios_base::app);

    for (config=0; config<n_configs; config++)
    {
        ndump << n_1vec[config] << "\n" << n0vec[config] << "\n" << n1vec[config] << "\n" << 
                 n2vec[config] << "\n" << n3vec[config] << "\n" <<
                 n4vec[config] << "\n" << n5vec[config] << "\n" <<
                 n6vec[config] << "\n";
    }

    ndump.close();

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

blocking::~blocking()
{
}
