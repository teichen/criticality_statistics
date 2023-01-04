// Gofr.h
#ifndef _GOFR
#define _GOFR

#include "SimSpace.h"

using namespace std;

class Gofr
{

public:

    bool mem_test;

    Gofr();

    SimSpace lattice; // define the coarse lattice
    int L;   // length of lattice (number of sites)
    int dim; // dimensionality of lattice
    int nL;

    void rdf(int*,double*);
    void initarrays();

    double dr,rL;
    int rbins;

    int* sites;

    double* r;
    int* nhist;
    double* n_ave;
    double* bin_vols;
    double* n_ideal;

    ~Gofr(); 

private:

};

#endif

