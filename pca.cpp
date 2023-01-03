#include "pca.h"
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <math.h>

#include <algorithm>

#include <sys/time.h>
#include <unistd.h>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <sstream>
using std::ifstream;
#include <cstring>

#include <omp.h>

#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

const double PI = 3.1415926535897932384626433832795028841971693;

using namespace std;

pca::pca()
{
    std::ostringstream filenameStream;
    std:string filename;

    ifstream fin;
    string line;

    std::ofstream tdump;
    std::ofstream edump;

    int i,j,k;
    int ii,jj,kk;

    int L[2];

    L[0] = 8;
    L[1] = 16;

    int n_configs;
    n_configs = 10;;

    int config;
    double nave;

    // npts sampled at equilibrium

    double n[n_configs*L[1]*L[1]*L[1]];
    double n_trans[n_configs*L[1]*L[1]*L[1]];

    int s;

    for (i=1; i<=n_configs; i++)
    {
        // read in data

        filenameStream << "./L" << L[1] << "/netlib/n" << i << ".dat";
        filename = filenameStream.str();

        fin.open(filename.c_str());

        k = (i-1)*L[1]*L[1]*L[1];
        while( getline(fin,line) )
        {
            n[k] = atoi (line.c_str());
            k++;
        }

        fin.close();
        filenameStream.str("");
    }

    // shift relative to site means

    for (j=0; j<(L[1]*L[1]*L[1]); j++)
    {
        nave = 0.0;
        for (i=1; i<=n_configs; i++)
        {
            nave = nave + n[(i-1)*L[1]*L[1]*L[1]+j];
        }

        nave = nave / (double)(n_configs);

        for (i=1; i<=n_configs; i++)
        {
            n[(i-1)*L[1]*L[1]*L[1]+j] = n[(i-1)*L[1]*L[1]*L[1]+j] - nave;
        }
    }

    // calculate transpose of data

    for (j=0; j<(L[1]*L[1]*L[1]); j++)
    {
        for (i=1; i<=n_configs; i++)
        {
            n_trans[j*n_configs+(i-1)] = n[(i-1)*L[1]*L[1]*L[1]+j];
        }
    }

    // PCA (SVD algorithm)

    // PCA of transpose of data

    gsl_matrix_view n_gsl = gsl_matrix_view_array (n_trans, L[1]*L[1]*L[1], n_configs);

    gsl_vector *s_svd = gsl_vector_alloc (n_configs);
    gsl_matrix *v_svd = gsl_matrix_alloc (n_configs, n_configs);
    gsl_vector *w_svd = gsl_vector_alloc (n_configs);

    gsl_linalg_SV_decomp (&n_gsl.matrix, v_svd, s_svd, w_svd);

    filenameStream << "./pcvals.dat";
    filename = filenameStream.str();
    edump.open(filename.c_str(), std::ios_base::app);

    for (j=0; j<n_configs; j++)
    {
        edump << gsl_vector_get(s_svd, j) << "\n";
    }

    edump.close();
    filenameStream.str("");

    // calculate principal components from transpose space
    // n = U * S * transpose(V)
    // transpose(n) = V * S * transpose(U)
    // V = transpose(n) * U * S
    // transpose(V) = S * transpose(U) * n

    filenameStream << "./pccomps.dat";
    filename = filenameStream.str();
    edump.open(filename.c_str(), std::ios_base::app);

    for (i=0; i<n_configs; i++)
    {
        for (j=0; j<L[1]; j++)
        {
            vt = 0.0;
            for (k=0; k<n_configs; k++)
            {
                vt = vt + gsl_vector_get(s_svd, i)*gsl_matrix_get(v_svd, k, i)*n[k*L[1]*L[1]*L[1]+j];
            }
            // print out transpose(V)
            edump << vt << "\n";
        }
    }

    edump.close();
    filenameStream.str("");

    gsl_vector_free (s_svd);
    gsl_matrix_free (v_svd);
    gsl_vector_free (w_svd);
}

pca::~pca()
{
}
