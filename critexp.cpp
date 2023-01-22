/*
    The critexp class calculates the critical exponents
    from the near-equilibrium statistics.
*/
#include "critexp.h"
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

critexp::critexp()
{
    int b;
    b = 4;

    L[0] = 8;
    L[1] = 16;

    n_configs = 10;

    jtypes = 5;
    jnum   = 8;

    mem_test = false;
    initarrays();

    int i,j,k;
    int ii,jj,kk;

    // read in data

    read_set_field(nb, b/2);
    read_set_field(na, b);

    // calculate averages
    double na_ave[jnum];
    double nb_ave[jnum];
    calc_averages(na, na_ave);
    calc_averages(nb, nb_ave);

    // calculate correlation functions
    double a[jnum*jtypes*jnum*jtypes];
    double across[jnum*jtypes*jnum*jtypes];

    calc_correlations(na, na, na_ave, na_ave, a);
    calc_correlations(na, nb, na_ave, nb_ave, across);

    write_correlations(a, "./a.dat");
    write_correlations(across, "./across.dat");

    // calculate condition numbers

    // calculate stability matrix
    // Ax=B --> x = linsolve(A,B)
    // t = linsolve(a,across)

    double t[jnum * jtypes * jnum * jtypes];
    double t_i[jnum * jtypes];
    double across_i[jnum * jtypes];
 
    for (i=0; i<(jnum*jtypes); i++)
    {
        for (j=0; j<(jnum*jtypes); j++)
        {
            across_i[j] = across[j*(jnum*jtypes)+i];
        }

        stability_matrix(a, across_i, t_i);

        write_t_row(t_i, "./t_stab.dat");

        for (j=0; j<(jnum*jtypes); j++)
        {
            t[j*(jnum*jtypes)+i] = t_i[j];
        }
    }

    // condition number:
    // int gsl_linalg_cholesky_rcond (const gsl_matrix * cholesky, double * rcond, gsl_vector * work)

    // eigensystem:
    double eval[jnum * jtypes];
    double evec[jnum * jtypes * jnum * jtypes];

    eigensystem(t, eval, evec);

    write_exponents(eval, "./exponent.dat");

}

void critexp::read_set_field(double* n, int b)
{    
    std::ostringstream filenameStream;
    std:string filename;

    ifstream fin;
    string line;

    int i,j;
    for (i=1; i<=n_configs; i++)
    {
        for (j=0; j<jnum; j++)
        {
            n[(i-1)*jnum+j] = 0.0;
        }

        filenameStream << "./L" << L[1] << "/" << i << "/nvec_b" << b << ".dat";
        filename = filenameStream.str();

        fin.open(filename.c_str());

        j = (i-1)*jnum;
        while( getline(fin,line) )
        {
            n[j] = atof (line.c_str());

            j++;
        }

        fin.close();
        filenameStream.str("");
    }
}

void critexp::write_correlations(double* a, string filename)
{
    std::ofstream edump;
    edump.open(filename, std::ios_base::app);

    int j,k;
    for (j=0; j<(jnum*jtypes); j++)
    {
        for (k=0; k<(jnum*jtypes); k++)
        {
            edump << a[j*jnum*jtypes+k] << "\n";
        }
    }

    edump.close();
}

void critexp::write_t_row(double* t_i, string filename)
{
    std::ofstream tdump;
    tdump.open(filename, std::ios_base::app);

    int j;
    for (j=0; j<(jnum*jtypes); j++)
    {
        tdump << t_i[j] << "\n";
    }

    tdump.close();
}

void critexp::calc_averages(double* n, double* nave)
{
    int i,j;
    for (i=0; i<jnum; i++)
    {
        nave[i] = 0.0;

        for (j=0; j<n_configs; j++)
        {
            nave[i] = nave[i] + n[j*jnum+i];
        }

        nave[i] = nave[i] / (double)(n_configs);
    }
}

void critexp::calc_correlations(double* n1, double* n2, double* n1_ave, double* n2_ave, double* c12)
{
    int i,j,k;

    for (i=0; i<jnum; i++)
    {
        for (j=0; j<jnum; j++)
        {
            c12[i*jnum+j] = 0.0;

            for (k=0; k<n_configs; k++)
            {
                c12[i*jnum+j] = c12[i*jnum+j] + na[k*jnum+i]*na[k*jnum+j];
            }
            c12[i*jnum+j] = c12[i*jnum+j] / (double)(n_configs);

            c12[i*jnum+j] = c12[i*jnum+j] - n1_ave[i]*n2_ave[j];
        }
    }
}

void critexp::stability_matrix(double* a, double* across_i, double* t)
{
    gsl_matrix_view a_gsl = gsl_matrix_view_array (a, jnum*jtypes, jnum*jtypes);
    gsl_vector_view across_gsl_i = gsl_vector_view_array (across_i, (jnum*jtypes));

    gsl_vector *t_i = gsl_vector_alloc (jnum*jtypes);

    gsl_permutation * p = gsl_permutation_alloc (jnum*jtypes);

    gsl_linalg_LU_decomp (&a_gsl.matrix, p, &s);
    gsl_linalg_LU_solve (&a_gsl.matrix, p, &across_gsl_i.vector, t_i);

    int j;
    for (j=0; j<(jnum*jtypes); j++)
    {
        t[j] = gsl_vector_get (t_i, j);
    }
    gsl_permutation_free (p);
    gsl_vector_free (t_i);
}

void critexp::eigensystem(double* t, double* eval, double* evec)
{
    gsl_matrix_view t_gsl = gsl_matrix_view_array (t, (jnum*jtypes), (jnum*jtypes));

    gsl_vector *s_svd = gsl_vector_alloc (jnum*jtypes);
    gsl_matrix *v_svd = gsl_matrix_alloc (jnum*jtypes, jnum*jtypes);
    gsl_vector *w_svd = gsl_vector_alloc (jnum*jtypes);

    gsl_vector_complex *gsl_eval = gsl_vector_complex_alloc (jnum*jtypes);
    gsl_matrix_complex *gsl_evec = gsl_matrix_complex_alloc (jnum*jtypes, jnum*jtypes);

    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (jnum*jtypes);

    gsl_eigen_nonsymmv (&t_gsl.matrix, gsl_eval, gsl_evec, w);

    gsl_linalg_SV_decomp (&t_gsl.matrix, v_svd, s_svd, w_svd);

    gsl_eigen_nonsymmv_free (w);

    gsl_eigen_nonsymmv_sort (gsl_eval, gsl_evec, GSL_EIGEN_SORT_ABS_ASC);

    gsl_complex eval_i;
    gsl_complex evec_ij;

    int i,j;
    for (i=0; i<jnum*jtypes; i++)
    {
        eval_i  = gsl_vector_complex_get (gsl_eval, i);
        eval[i] = (double)(GSL_REAL(eval_i));

        for (j=0; j<jnum*jtypes; j++)
        {
            evec_ij = gsl_matrix_complex_get (gsl_evec, i, j);
            evec[i*jnum*jtypes + j] = (double)(GSL_REAL(evec_ij));
        }
    }

    gsl_vector_free (s_svd);
    gsl_matrix_free (v_svd);
    gsl_vector_free (w_svd);

    gsl_vector_complex_free (gsl_eval);
    gsl_matrix_complex_free (gsl_evec);
}

void critexp::write_exponents(double* eval, string filename)
{
    std::ofstream edump;
    edump.open(filename, std::ios_base::app);

    int i;
    for (i=0; i<(jnum*jtypes); i++)
    {
        edump.open(filename, std::ios_base::app);

        if (eval[i]>0.0)
        {
            edump << log(eval[i])/log(2) << "\n";
        }

        edump.close();
    }
}

void critexp::initarrays()
{
    na       = (double*) calloc (n_configs * jnum, sizeof(double));
    nb       = (double*) calloc (n_configs * jnum, sizeof(double));

    mem_test = true;
}

critexp::~critexp()
{
    if(mem_test==true)
    {
    delete [] na;
    delete [] nb;

    cout << "Deallocate critexp memory" << endl;

    }
}
