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
    int L; // length of lattice (number of sites)
    int Lfine; // fine grid
    int dim; // dimensionality of lattice

    int nL;
    int nLfine;

    int bshift;

    int nvols;

    int get_cell(int,int,int,int,int,int,int);
    double separation(double[3],double[3]);

    ~blocking(); 

private:

};

#endif

