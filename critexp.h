// critexp.h
#ifndef _CRITEXP
#define _CRITEXP

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

    void calc_averages(double*, double*);
    void calc_correlations(double*, double*, double*, double*, double*);

    double* t;
    double* across_i;
    int s;

    void initarrays();

    void read_set_field(int*, int);

    ~critexp(); 

private:

};

#endif

