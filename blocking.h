// blocking.h
#ifndef _BLOCKING
#define _BLOCKING

#include "SimSpace.h"

using namespace std;

class blocking
{

public:

    bool mem_test;

    blocking();
    void coarse_grain_field(int*, int);

    SimSpace lattice; // define the coarse lattice
    int L;   // length of lattice (number of sites)
    int dim; // dimensionality of lattice

    SimSpace coarse_lattice;
    int Lb;

    int nL;

    int b;
    int bshift;

    double n_collective[8];

    double* nb;

    int nsites;

    void lattice_map(int, int*);

    void initarrays();

    ~blocking(); 

private:

};

#endif

