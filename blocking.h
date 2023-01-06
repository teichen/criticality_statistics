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
    blocking(int);

    SimSpace lattice; // define the coarse lattice
    int L;   // length of lattice (number of sites)
    int dim; // dimensionality of lattice

    SimSpace coarse_lattice;
    int Lb;

    int nL;

    int bshift;

    int nsites;

    int get_cell(int,int,int,int,int,int,int);

    ~blocking(); 

private:

};

#endif

