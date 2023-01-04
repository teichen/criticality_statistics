#include <iostream>
#include <cassert>
using std::cerr;
using std::cout;
using std::endl;
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "../Gofr.h"

using namespace std;
int main()
{
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
        n[i] = 2*(rand() % 2)-1;
    }

    gofr.rdf(n, g);

    for (i=0; i<gofr.rbins; i++)
    {
        cout << g[i] << endl;
    }
}

