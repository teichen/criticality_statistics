// SimSpace.h
#ifndef _SIMSPACE
#define _SIMSPACE

#include <iostream>

using namespace std;

class SimSpace
{

public:

    bool mem_test;

    SimSpace(); 

    int L; // length of lattice (number of sites)
    int dim; // dimensionality of lattice

    double separation(double[3],double[3]);

    int* unpack_position(int);

    int flatten_position(int, int, int);

    int* nearest_neighbors(int*);
    int* nearest_neighbor_values(int*, int*);

    ~SimSpace(); 

private:

};

#endif
