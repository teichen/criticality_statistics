// critexp.h
#ifndef _CRITEXP
#define _CRITEXP

#include <gsl/gsl_math.h>

#include <iostream>

using namespace std;

class critexp
{

public:

    bool mem_test;

    critexp();

    int L[2];
    int n_configs;
    int jtypes,jnum;

    double* na;
    double* nb;

    void read_set_field(double*, int);
    void write_correlations(double*, string);
    void write_t_row(double*, string);
    void calc_averages(double*, double*);
    void calc_correlations(double*, double*, double*, double*, double*);
    void stability_matrix(double*, double*, double*);
    void eigensystem(double*, double*, double*);
    void write_exponents(double*, string);

    int s;

    void initarrays();

    void read_set_field(int*, int);

    ~critexp(); 

private:

};

#endif

