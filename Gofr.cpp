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
    L = lattice.L; // length of lattice (number of sites)
    Lfine = lattice.Lfine; // fine grid

    dim = lattice.dim; // dimensionality of lattice

    nL = pow(L,dim);

    rbins = 100+1; // resolution of rdf

    mem_test = false;
    initarrays();

}

void Gofr::vaprdf(int* n, double* gofr_vap)
{
    int i,j;

    double ri[3];
    double rj[3];
    double d[3];
    double rmag;
    int bin_ij;
	
    double volume,density;

    int dimR;

    dr = L/(double)(rbins-1);

    for (i=0; i<rbins; i++)
    {
        r[i] = i*dr;
        nhist[i] = 0;
        n_ave[i] = 0.0;
        bin_vols[i] = 0.0;
        n_ideal[i] = 0.0;
    }

    int nvap;
    nvap = 0;

    for (i=0; i<nL; i++)
    {
        if (n[i]==0)
        {
            nvap++;
        }
    }

    // create an array of the vapor sites

    int* vapsites;
    vapsites = (int*) calloc (nvap, sizeof(int));

    j = 0;
    for (i=0; i<nL; i++)
    {
        if (n[i]==0)
        {
            vapsites[j] = i;
            j++;
        }
    }

    // calculate radial distribution function for the vapor sites
    for (j=0; j<rbins; j++)
    {
        gofr_vap[j] = 0.0;
    }

    int ii,jj;

    for (i=0; i<(int)(nvap-1); i++)
    {
        ii = vapsites[i];

        ri[2] = (int)(ii % L) + 0.5;
        ri[1] = ((int)((int)(ii-(int)(ii % L))/L) % L) + 0.5;
        ri[0] = ((int)((int)(ii-(int)(ii % L)-(int)((int)(ii-(int)(ii % L)/L) % L)*L)/pow(L,2)) % L) + 0.5;

        for (j=(i+1); j<nvap; j++)
        {
            jj = vapsites[j];

            rj[2] = (int)(jj % L) + 0.5;
            rj[1] = ((int)((int)(jj-(int)(jj % L))/L) % L) + 0.5;
            rj[0] = ((int)((int)(jj-(int)(jj % L)-(int)((int)(jj-(int)(jj % L)/L) % L)*L)/pow(L,2)) % L) + 0.5;

            d[0] = 0.0;
            d[1] = 0.0;
            d[2] = 0.0;

            for (dimR=0; dimR<dim; dimR++)
            {
                if ((ri[dimR]-rj[dimR])>(double)(0.5*L))
                {
                    d[dimR] = ri[dimR] - rj[dimR] - L;
                }
                else if((ri[dimR]-rj[dimR])<(double)(-0.5*L))
                {
                    d[dimR] = ri[dimR] - rj[dimR] + L;
                }
                else
                {
                    d[dimR] = ri[dimR] - rj[dimR];
                }
            }
            rmag = sqrt(pow(d[0],2) + pow(d[1],2) + pow(d[2],2));

            bin_ij = ceil(rmag/dr);

            if ( (bin_ij <= rbins) && (bin_ij > 0) )
            {
                nhist[bin_ij-1] = nhist[bin_ij-1] + 1;
            }
        }
    }

    volume = nL;
    density = nvap/(double)volume;

    for (i=0; i<rbins; i++)
    {
        n_ave[i] = 2*nhist[i]/(double)nvap;
        bin_vols[i] = (4*PI/3)*( pow(r[i]+dr,3) - pow(r[i],3) );
        n_ideal[i] = density*bin_vols[i];
        gofr_vap[i] = n_ave[i]/n_ideal[i];

    }

    delete [] vapsites;

}

void Gofr::liqrdf(double* v, double* gofr_liq)
{
    int i,j;

    double ri[3];
    double rj[3];
    double d[3];
    double rmag;
    int bin_ij;
	
    double volume,density;

    int dimR;

    dr = L/(double)(rbins-1);

    for (i=0; i<rbins; i++)
    {
        r[i] = i*dr;
        nhist[i] = 0;
        n_ave[i] = 0.0;
        bin_vols[i] = 0.0;
        n_ideal[i] = 0.0;
    }

    int nliq;
    nliq = 0;

    for (i=0; i<nL; i++)
    {
        if (v[i]<=0.25)
        {
            nliq++;
        }
    }

    // create an array of the liquid sites

    int* liqsites;
    liqsites = (int*) calloc (nliq, sizeof(int));

    j = 0;
    for (i=0; i<nL; i++)
    {
        if (v[i]<=0.25)
        {
            liqsites[j] = i;
            j++;
        }
    }

    // calculate radial distribution function for the liquid sites
    for (j=0; j<rbins; j++)
    {
        gofr_liq[j] = 0.0;
    }

    int ii,jj;

    for (i=0; i<(int)(nliq-1); i++)
    {
        ii = liqsites[i];

        ri[2] = (int)(ii % L) + 0.5;
        ri[1] = ((int)((int)(ii-(int)(ii % L))/L) % L) + 0.5;
        ri[0] = ((int)((int)(ii-(int)(ii % L)-(int)((int)(ii-(int)(ii % L)/L) % L)*L)/pow(L,2)) % L) + 0.5;

        for (j=(i+1); j<nliq; j++)
        {
            jj = liqsites[j];

            rj[2] = (int)(jj % L) + 0.5;
            rj[1] = ((int)((int)(jj-(int)(jj % L))/L) % L) + 0.5;
            rj[0] = ((int)((int)(jj-(int)(jj % L)-(int)((int)(jj-(int)(jj % L)/L) % L)*L)/pow(L,2)) % L) + 0.5;

            d[0] = 0.0;
            d[1] = 0.0;
            d[2] = 0.0;

            for (dimR=0; dimR<dim; dimR++)
            {
                if ((ri[dimR]-rj[dimR])>(double)(0.5*L))
                {
                    d[dimR] = ri[dimR] - rj[dimR] - L;
                }
                else if((ri[dimR]-rj[dimR])<(double)(-0.5*L))
                {
                    d[dimR] = ri[dimR] - rj[dimR] + L;
                }
                else
                {
                    d[dimR] = ri[dimR] - rj[dimR];
                }
            }
            rmag = sqrt(pow(d[0],2) + pow(d[1],2) + pow(d[2],2));

            bin_ij = ceil(rmag/dr);

            if ( (bin_ij <= rbins) && (bin_ij > 0) )
            {
                nhist[bin_ij-1] = nhist[bin_ij-1] + 1;
            }
        }
    }

    volume = nL;
    density = nliq/(double)volume;

    for (i=0; i<rbins; i++)
    {
        n_ave[i] = 2*nhist[i]/(double)nliq;
        bin_vols[i] = (4*PI/3)*( pow(r[i]+dr,3) - pow(r[i],3) );
        n_ideal[i] = density*bin_vols[i];
        gofr_liq[i] = n_ave[i]/n_ideal[i];

    }

    delete [] liqsites;

}

void Gofr::initarrays()
{
    r = (double*) calloc (rbins, sizeof(double));
    nhist = (int*) calloc (rbins, sizeof(int));
    n_ave = (double*) calloc (rbins, sizeof(double));
    bin_vols = (double*) calloc (rbins, sizeof(double));
    n_ideal = (double*) calloc (rbins, sizeof(double));

    mem_test = true;
}

Gofr::~Gofr()
{
    if(mem_test=true)
    {
    delete [] r;
    delete [] nhist;
    delete [] n_ave;
    delete [] bin_vols;
    delete [] n_ideal;
    }
}
