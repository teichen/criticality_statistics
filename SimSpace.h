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

    void set_dimensions(int);

    double separation(int*, int*);

    void unpack_position(int, int[3]);

    int flatten_position(int, int, int);

    void nearest_neighbors(int*, int[3]);
    void diagonal_in_plane(int*, int[6]);
    void cubic_diagonal(int*, int[4]);
    void principal_planes(int*, int[9]);
    void diagonal_planes(int*, int[18]);
    void tetrahedral_vertices(int*, int[6]);
    void next_nearest_neighbors(int*, int[6]);

    ~SimSpace(); 

private:

};

#endif
