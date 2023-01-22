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

    // critexp.cpp unit tests
    
    // pca.cpp unit tests

}

