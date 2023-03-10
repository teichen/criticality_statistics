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
#include <gsl/gsl_vector.h>

const double PI = 3.1415926535897932384626433832795028841971693;

using namespace std;

pca::pca()
{
    /* principal component analysis
    */
    L[0] = 8;
    L[1] = 16;

    n_configs = 10;

    mem_test = false;
    initarrays();

    int i,j,k;
    int ii,jj,kk;

    int config;

    // read in data
    read_set_field(n);

    // shift relative to site means
    mean_shift_field(n);

    // calculate transpose of data
    transpose_field(n, ntrans);

    // PCA (SVD algorithm)
    // PCA of transpose of data

    // principal values
    double s[n_configs];
    // transpose of V
    double vt[n_configs * L[1]];

    principal_component_analysis(n, ntrans, s, vt);

    write_principal_values(s, "./pcvals.dat");
}

void pca::read_set_field(double* n)
{    
    std::ostringstream filenameStream;
    std:string filename;

    ifstream fin;
    string line;

    int i,j;
    for (i=1; i<=n_configs; i++)
    {
        filenameStream << "./L" << L[1] << "/netlib/n" << i << ".dat";
        filename = filenameStream.str();

        fin.open(filename.c_str());

        j = (i-1)*L[1]*L[1]*L[1];
        while( getline(fin,line) )
        {
            n[j] = atoi (line.c_str());
            j++;
        }

        fin.close();
        filenameStream.str("");
    }
}

void pca::mean_shift_field(double* n)
{
    double nave;
    int i,j;

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
}

void pca::transpose_field(double* n, double* ntrans)
{
    int i,j;
    for (j=0; j<(L[1]*L[1]*L[1]); j++)
    {
        for (i=1; i<=n_configs; i++)
        {
            ntrans[j*n_configs+(i-1)] = n[(i-1)*L[1]*L[1]*L[1]+j];
        }
    }
}

void pca::principal_component_analysis(double* n, double* ntrans, double* s, double* vt)
{
    /* principal component analysis
       Args:
               n (double*)     :
               ntrans (double*):
               s (double*)     :
               vt (double*)    :
    */
    gsl_matrix_view n_gsl = gsl_matrix_view_array (ntrans, L[1]*L[1]*L[1], n_configs);

    gsl_vector *s_svd = gsl_vector_alloc (n_configs);
    gsl_matrix *v_svd = gsl_matrix_alloc (n_configs, n_configs);
    gsl_vector *w_svd = gsl_vector_alloc (n_configs);

    gsl_linalg_SV_decomp (&n_gsl.matrix, v_svd, s_svd, w_svd);

    int i,j,k;
    for (i=0; i<n_configs; i++)
    {
        s[i] = gsl_vector_get(s_svd, i);
    }
    // calculate principal components from transpose space
    // n = U * S * transpose(V)
    // transpose(n) = V * S * transpose(U)
    // V = transpose(n) * U * S
    // transpose(V) = S * transpose(U) * n

    for (i=0; i<n_configs; i++)
    {
        for (j=0; j<L[1]; j++)
        {
            vt[i*L[1] + j] = 0.0;
            for (k=0; k<n_configs; k++)
            {
                vt[i*L[1] + j] += gsl_vector_get(s_svd, i) * gsl_matrix_get(v_svd, k, i) * n[k*L[1]*L[1]*L[1]+j];
            }
        }
    }
    gsl_vector_free (s_svd);
    gsl_matrix_free (v_svd);
    gsl_vector_free (w_svd);
}

void pca::write_principal_values(double* s, string filename)
{
    std::ofstream edump;
    edump.open(filename, std::ios_base::app);

    int j;
    for (j=0; j<n_configs; j++)
    {
        edump << s[j] << "\n";
    }

    edump.close();
}

void pca::initarrays()
{
    n       = (double*) calloc (n_configs * pow(L[1], 3), sizeof(double));
    ntrans  = (double*) calloc (n_configs * pow(L[1], 3), sizeof(double));

    mem_test = true;
}

pca::~pca()
{
    if(mem_test==true)
    {
    delete [] n;
    delete [] ntrans;

    cout << "Deallocate pca memory" << endl;

    }
}

