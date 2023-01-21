// pca.h
#ifndef _PCA
#define _PCA

#include <iostream>

using namespace std;

class pca
{

public:

    bool mem_test;

    pca();

    int L[2];
    int n_configs;

    double* n;
    double* ntrans;

    void initarrays();

    void read_set_field(double*);
    void mean_shift_field(double*);
    void transpose_field(double*, double*);
    void principal_component_analysis(double*, double*, double*, double*);
    void write_principal_values(double*, string);

    ~pca(); 

private:

};

#endif

