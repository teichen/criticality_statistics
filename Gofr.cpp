/*
    The Gofr class calculates the radial distribution functions
    of a binary discretized field.
*/
#include "Gofr.h"
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

#include <time.h>
#include <math.h>

const double PI = 3.1415926535897932384626433832795028841971693;

using namespace std;

Gofr::Gofr()
{
    L   = lattice.L;   // length of lattice (number of sites)
    dim = lattice.dim; // dimensionality of lattice
    nL  = pow(L,dim);

    rbins = 100+1; // resolution of rdf

    mem_test = false;
    initarrays();
}

void Gofr::rdf(int* n, double* gofr)
{
    int i,j;

    int* ri;
    int* rj;

    double separation;
    int bin_ij;
	
    double volume, density;

    int dimR;

    dr = L / (double)(rbins - 1);

    for (i=0; i<rbins; i++)
    {
        r[i]     = i*dr;
        nhist[i] = 0;
        n_ave[i] = 0.0;
        bin_vols[i] = 0.0;
        n_ideal[i]  = 0.0;
    }

    sites = (int*) calloc (nL, sizeof(int));

    int nsites;

    j = 0;
    for (i=0; i<nL; i++)
    {
        if (n[i]==0)
        {
            sites[j] = i;
            j++;
        }
    }
    nsites = j;

    // calculate radial distribution function
    for (i=0; i<rbins; i++)
    {
        gofr[i] = 0.0;
    }

    int ii,jj;

    for (i=0; i<(int)(nsites - 1); i++)
    {
        ii = sites[i];

        ri = lattice.unpack_position(ii);

        for (j=(i+1); j<(int)(nsites - 1); j++)
        {
            jj = sites[j];

            rj = lattice.unpack_position(jj);

            separation = lattice.separation(ri, rj);

            bin_ij = ceil(separation / dr);

            if ( (bin_ij <= rbins) && (bin_ij > 0) )
            {
                nhist[bin_ij-1] = nhist[bin_ij-1] + 1;
            }
        }
    }

    volume = nL;
    density = nsites / (double)volume;

    for (i=0; i<rbins; i++)
    {
        n_ave[i]    = 2 * nhist[i] / (double)nsites;
        bin_vols[i] = (4*PI/3) * ( pow(r[i] + dr, 3) - pow(r[i], 3) );
        n_ideal[i]  = density * bin_vols[i];
        gofr[i]     = n_ave[i] / n_ideal[i];

    }
    delete [] sites;
}

void Gofr::initarrays()
{
    r     = (double*) calloc (rbins, sizeof(double));
    nhist = (int*) calloc (rbins, sizeof(int));
    n_ave = (double*) calloc (rbins, sizeof(double));
    bin_vols = (double*) calloc (rbins, sizeof(double));
    n_ideal  = (double*) calloc (rbins, sizeof(double));

    mem_test = true;
}

Gofr::~Gofr()
{
    if(mem_test==true)
    {
    delete [] r;
    delete [] nhist;
    delete [] n_ave;
    delete [] bin_vols;
    delete [] n_ideal;
    }
}
