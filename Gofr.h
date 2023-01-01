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
    int L; // length of lattice (number of sites)
    int Lfine; // fine grid
    int dim; // dimensionality of lattice
    int nL;

    void vaprdf(int*,double*);
    void liqrdf(double*,double*);
    void initarrays();

    double dr,rL;
    int rbins;
   
    double* r;
    int* nhist;
    double* n_ave;
    double* bin_vols;
    double* n_ideal;

    ~Gofr(); 

private:

};

#endif

