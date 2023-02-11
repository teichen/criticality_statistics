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
}

critexp::critexp(bool logger, int& configs, int& jnum, double* n1, double* n2)
{
    /* calculate critical exponents for the collective statistics of a number
       of configurations sampled near equilibrium but coarsened at different
       blocking levels, n1 and n2

       Args:
               logger (bool)      : if true, print Glauber dynamics
               configs (int)      : number of configurations
               jnum (int)         : number of collective statistics
               n1 (array)         : first array of collective statistics
               n2 (array)         : second array of collective statistics
    */
    logging   = logger;
    n_configs = configs;
    n_coord   = jnum;

    mem_test = false;
    initarrays();

    int i;
    for (i=0; i<n_coord; i++)
    {
        na[i] = n1[i];
        nb[i] = n2[i];
    }

    // TODO: read_set_field() from disk (collective coordinates averaged 
    //       with various blocking levels, b)
    // read_set_field(nb, b/2);
    // read_set_field(na, b);

    // calculate averages
    double na_ave[n_coord];
    double nb_ave[n_coord];
    calc_averages(na, na_ave);
    calc_averages(nb, nb_ave);

    // calculate correlation functions
    double a[n_coord * n_coord];
    double across[n_coord * n_coord];

    calc_correlations(na, na, na_ave, na_ave, a);
    calc_correlations(na, nb, na_ave, nb_ave, across);

    if(logging)
    {
        write_correlations(a, "./a.dat");
        write_correlations(across, "./across.dat");
    }

    /* calculate condition numbers

       calculate stability matrix
       Ax=B --> x = linsolve(A,B)
       t = linsolve(a,across)
    */
    double t[n_coord * n_coord];
    double t_i[n_coord];
    double across_i[n_coord];

    int j;
    for (i=0; i<n_coord; i++)
    {
        for (j=0; j<n_coord; j++)
        {
            across_i[j] = across[j*n_coord + i];
        }

        stability_matrix(a, across_i, t_i);

        write_t_row(t_i, "./t_stab.dat");

        for (j=0; j<n_coord; j++)
        {
            t[j*n_coord + i] = t_i[j];
        }
    }

    /* condition number:
       int gsl_linalg_cholesky_rcond (const gsl_matrix * cholesky, double * rcond, gsl_vector * work)
       eigensystem:
    */
    double eval[n_coord];
    double evec[n_coord * n_coord];

    eigensystem(t, eval, evec);

    for (i=0; i<n_coord; i++)
    {
        if (eval[i]>0.0)
        {
            exponents[i] = log(eval[i])/log(2);
        }
    }

    if(logging)
    {
        write_exponents(exponents, "./exponent.dat");
    }
}

void critexp::write_correlations(double* a, string filename)
{
    /* write out correlations
    */
    std::ofstream edump;
    edump.open(filename, std::ios_base::app);

    int j,k;
    for (j=0; j<n_coord; j++)
    {
        for (k=0; k<n_coord; k++)
        {
            edump << a[j*n_coord + k] << "\n";
        }
    }

    edump.close();
}

void critexp::write_t_row(double* t_i, string filename)
{
    /* write out the stability matrix
    */
    std::ofstream tdump;
    tdump.open(filename, std::ios_base::app);

    int j;
    for (j=0; j<n_coord; j++)
    {
        tdump << t_i[j] << "\n";
    }

    tdump.close();
}

void critexp::calc_averages(double* n, double* nave)
{
    /* calculate averages of the collective coordinates
    */
    int i,j;
    for (i=0; i<n_coord; i++)
    {
        nave[i] = 0.0;

        for (j=0; j<n_configs; j++)
        {
            nave[i] = nave[i] + n[j*n_coord+i];
        }

        nave[i] = nave[i] / (double)(n_configs);
    }
}

void critexp::calc_correlations(double* n1, double* n2, double* n1_ave, double* n2_ave, double* c12)
{
    /* calculate correlation between the two blocking levels
    */
    int i,j,k;

    for (i=0; i<n_coord; i++)
    {
        for (j=0; j<n_coord; j++)
        {
            c12[i*n_coord+j] = 0.0;

            for (k=0; k<n_configs; k++)
            {
                c12[i*n_coord+j] = c12[i*n_coord+j] + n1[k*n_coord+i]*n2[k*n_coord+j];
            }
            c12[i*n_coord+j] = c12[i*n_coord+j] / (double)(n_configs);

            c12[i*n_coord+j] = c12[i*n_coord+j] - n1_ave[i]*n2_ave[j];
        }
    }
}

void critexp::stability_matrix(double* a, double* across_i, double* t)
{
    /* calculate the stability matrix
    */
    gsl_matrix_view a_gsl = gsl_matrix_view_array (a, n_coord, n_coord);
    gsl_vector_view across_gsl_i = gsl_vector_view_array (across_i, (n_coord));

    gsl_vector *t_i = gsl_vector_alloc (n_coord);

    gsl_permutation * p = gsl_permutation_alloc (n_coord);

    gsl_linalg_LU_decomp (&a_gsl.matrix, p, &s);
    
    gsl_linalg_LU_solve (&a_gsl.matrix, p, &across_gsl_i.vector, t_i);

    int j;
    for (j=0; j<(n_coord); j++)
    {
        t[j] = gsl_vector_get (t_i, j);
    }
    gsl_permutation_free (p);
    gsl_vector_free (t_i);
}

void critexp::eigensystem(double* t, double* eval, double* evec)
{
    /* decompose the stability matrix into an eigensystem for the 
       calculation of critical exponents
    */
    gsl_matrix_view t_gsl = gsl_matrix_view_array (t, (n_coord), (n_coord));

    gsl_vector *s_svd = gsl_vector_alloc (n_coord);
    gsl_matrix *v_svd = gsl_matrix_alloc (n_coord, n_coord);
    gsl_vector *w_svd = gsl_vector_alloc (n_coord);

    gsl_vector_complex *gsl_eval = gsl_vector_complex_alloc (n_coord);
    gsl_matrix_complex *gsl_evec = gsl_matrix_complex_alloc (n_coord, n_coord);

    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (n_coord);

    gsl_eigen_nonsymmv (&t_gsl.matrix, gsl_eval, gsl_evec, w);

    gsl_linalg_SV_decomp (&t_gsl.matrix, v_svd, s_svd, w_svd);

    gsl_eigen_nonsymmv_free (w);

    gsl_eigen_nonsymmv_sort (gsl_eval, gsl_evec, GSL_EIGEN_SORT_ABS_ASC);

    gsl_complex eval_i;
    gsl_complex evec_ij;

    int i,j;
    for (i=0; i<n_coord; i++)
    {
        eval_i  = gsl_vector_complex_get (gsl_eval, i);
        eval[i] = (double)(GSL_REAL(eval_i));

        for (j=0; j<n_coord; j++)
        {
            evec_ij = gsl_matrix_complex_get (gsl_evec, i, j);
            evec[i*n_coord + j] = (double)(GSL_REAL(evec_ij));
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
    /* write out critical exponents
    */
    std::ofstream edump;
    edump.open(filename, std::ios_base::app);

    int i;
    for (i=0; i<(n_coord); i++)
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
    exponents = (double*) calloc (n_coord, sizeof(double));
    na        = (double*) calloc (n_coord, sizeof(double));
    nb        = (double*) calloc (n_coord, sizeof(double));

    mem_test  = true;
}

critexp::~critexp()
{
    if(mem_test==true)
    {
    delete [] exponents;
    delete [] na;
    delete [] nb;

    cout << "Deallocate critexp memory" << endl;

    }
}
