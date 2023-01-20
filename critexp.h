// critexp.h
#ifndef _CRITEXP
#define _CRITEXP

#include <gsl/gsl_math.h>

using namespace std;

class critexp
{

public:

    bool mem_test;

    critexp();

    int n_configs;
    int jtypes,jnum;

    int* na;
    int* nb;

    void write_correlations(double*, string);
    void write_t_row(double*, string);
    void calc_averages(double*, double*);
    void calc_correlations(double*, double*, double*, double*, double*);
    void stability_matrix(double*, double*, double*);
    void eigensystem(double*, gsl_vector_complex*, gsl_matrix_complex*);
    void write_exponents(gsl_vector_complex*, string);

    int s;

    void initarrays();

    void read_set_field(int*, int);

    ~critexp(); 

private:

};

#endif

