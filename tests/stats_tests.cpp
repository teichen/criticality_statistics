#include <iostream>
#include <cassert>
using std::cerr;
using std::cout;
using std::endl;
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "../Gofr.h"
#include "../blocking.h"
#include "../critexp.h"
#include "../pca.h"

using namespace std;
int main()
{
    // Gofr.cpp unit tests
    Gofr gofr;

    int n[gofr.nL];
    double g[gofr.rbins];

    srand (time(NULL)); // re-initialize random seed

    int i;
    
    for (i=0; i<gofr.rbins; i++)
    {
        g[i] = 0.0;
    }
    // example field: two spin up
    for (i=0; i<gofr.nL; i++)
    {
        //n[i] = 2*(rand() % 2)-1;
        n[i] = -1;
    }
    n[0] = 1;
    n[2] = 1;

    gofr.rdf(n, g);

    double dr;
    dr = gofr.L / (double)(gofr.rbins - 1);

    int bin_02;
    bin_02 = ceil(2.0 / dr);

    for (i=0; i<gofr.rbins; i++)
    {
        if(i==(int)(bin_02-1))
        {
            assert(g[i]>0);
        }
        else
        {
            assert(g[i]==0);
        }
    }

    // blocking.cpp unit tests
    int b = 4; 

    blocking cg_lattice_b4;
    cg_lattice_b4.coarse_grain_field(n, b);

    int nLb;
    nLb = pow(cg_lattice_b4.Lb, 3);
    
    // validate the sum, example field has two spin up
    double sum_b = 0.0;
    sum_b = -1 * (pow(b, 3)-2) + 2;
    sum_b /= pow(b, 3);
    for (i=1; i<nLb; i++)
    {
        sum_b += -1;
    }
    
    assert(cg_lattice_b4.n_collective[0]==sum_b);

    // critexp.cpp unit tests
    bool logging;
    logging = 0;

    int n_configs = 2;
    int jnum      = 8;
    double n_b2[n_configs * jnum];
    double n_b4[n_configs * jnum];

    b = 2;
    blocking cg_lattice_b2;
    cg_lattice_b2.coarse_grain_field(n, b);
    for (i=0; i<jnum; i++)
    {
        n_b2[i] = cg_lattice_b2.n_collective[i];
    }
    b = 4;
    cg_lattice_b4.coarse_grain_field(n, b);
    for (i=0; i<jnum; i++)
    {
        n_b4[i] = cg_lattice_b4.n_collective[i];
    }

    // make some changes to the input field to generate new statistics
    n[0] = -1;
    n[2] = 1;
    n[4] = 1;

    b = 2;
    cg_lattice_b2.coarse_grain_field(n, b);
    for (i=0; i<jnum; i++)
    {
        n_b2[i+jnum] = cg_lattice_b2.n_collective[i];
    }
    b = 4;
    cg_lattice_b4.coarse_grain_field(n, b);
    for (i=0; i<jnum; i++)
    {
        n_b4[i+jnum] = cg_lattice_b4.n_collective[i];
    }

    critexp criticality(logging, n_configs, jnum, n_b2, n_b4);

    for (i=0; i<jnum; i++)
    {
        cout << criticality.exponents[i] << endl;
    }

    // TODO add assert();

    // pca.cpp unit tests

}

