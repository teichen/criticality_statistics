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

const double PI = 3.1415926535897932384626433832795028841971693;

using namespace std;

critexp::critexp()
{
    int b;

    b = 4;

    int jtypes,jnum;

    jtypes = 5;
    jnum = 8;

    double jvec[jtypes*jnum];

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
    n_configs = 10;

    // n_configs sampled at equilibrium
    double na[n_configs*jnum];
    double nb[n_configs*jnum];

    // initialize jnum*2 matrices for nave averages
    double nave[jnum*2];

    // initialize jnum*jnum matrices for nk,nj correlation functions
    double n_n[jnum*jnum];

    double a[jnum*jtypes*jnum*jtypes];
    double across[jnum*jtypes*jnum*jtypes];

    double t[jnum*jtypes*jnum*jtypes];

    double across_i[jnum*jtypes];
    int s;

    // read in data

    for (i=0; i<(jtypes*jnum); i++)
    {
        jvec[i] = 0.0;
    }

    jind = 0;

    for (i=1; i<=n_configs; i++)
    {
        for (j=0; j<jnum; j++)
        {
            na[(i-1)*jnum+j] = 0.0;
        }
        for (j=0; j<jnum; j++)
        {
            nb[(i-1)*jnum+j] = 0.0;
        }

        filenameStream << "./L" << L[1] << "/" << i << "/nvec_b" << b/2 << ".dat";
        filename = filenameStream.str();

        fin.open(filename.c_str());

        j = (i-1)*jnum;
        while( getline(fin,line) )
        {
            nb[j] = atof (line.c_str());

            j++;
        }

        fin.close();
        filenameStream.str("");

        filenameStream << "./L" << L[1] << "/" << i << "/nvec_b" << b << ".dat";
        filename = filenameStream.str();

        fin.open(filename.c_str());

        j = (i-1)*jnum;
        while( getline(fin,line) )
        {
            na[j] = atof (line.c_str());

            j++;
        }

        fin.close();
        filenameStream.str("");
    }

    // calculate averages
    for (i=0; i<jnum; i++)
    {
        nave[i*2+0] = 0.0;
        nave[i*2+1] = 0.0;

        for (j=0; j<(n_configs*npts[1]); j++)
        {
            nave[i*2+0] = nave[i*2+0] + na[j*jnum+i];
        }
        for (j=0; j<(n_configs*npts[1]); j++)
        {
            nave[i*2+1] = nave[i*2+1] + nb[j*jnum+i];
        }

        nave[i*2+0] = nave[i*2+0] / (double)(n_configs*npts[1]);
        nave[i*2+1] = nave[i*2+1] / (double)(n_configs*npts[1]);
    }

    // calculate correlation functions
    for (i=0; i<jnum; i++)
    {
        for (j=0; j<jnum; j++)
        {
            n_n[i*jnum+j] = 0.0;

            for (k=0; k<n_configs; k++)
            {
                n_n[i*jnum+j] = n_n[i*jnum+j] + na[k*jnum+i]*na[k*jnum+j];
            }
            n_n[i*jnum+j] = n_n[i*jnum+j] / (double)(n_configs);

            n_n[i*jnum+j] = n_n[i*jnum+j] - nave[i*2+1]*nave[j*2+1];
        }
    }

    for (i=0; i<jnum; i++)
    {
        for (j=0; j<jnum; j++)
        {
            a[i*jnum*jtypes+j] = n_n[i*jnum+j];
        }
    }

    for (i=0; i<jnum; i++)
    {
        for (j=0; j<jnum; j++)
        {
            n_n[i*jnum+j] = 0.0;

            for (k=0; k<n_configs; k++)
            {
                n_n[i*jnum+j] = n_n[i*jnum+j] + na[k*jnum+i]*nb[k*jnum+j];
            }
            n_n[i*jnum+j] = n_n[i*jnum+j]/(double)(n_configs);

            n_n[i*jnum+j] = n_n[i*jnum+j] - nave[i*2+1]*nave[j*2+0];
        }
    }

    for (i=0; i<jnum; i++)
    {
        for (j=0; j<jnum; j++)
        {
            across[i*jnum*jtypes+j] = n_n[i*jnum+j];
        }
    }

    filenameStream << "./a.dat";
    filename = filenameStream.str();
    edump.open(filename.c_str(), std::ios_base::app);

    for (j=0; j<(jnum*jtypes); j++)
    {
        for (k=0; k<(jnum*jtypes); k++)
        {
            edump << a[j*jnum*jtypes+k] << "\n";
        }
    }

    edump.close();
    filenameStream.str("");

    filenameStream << "./across.dat";
    filename = filenameStream.str();
    edump.open(filename.c_str(), std::ios_base::app);

    for (j=0; j<(jnum*jtypes); j++)
    {
        for (k=0; k<(jnum*jtypes); k++)
        {
            edump << across[j*jnum*jtypes+k] << "\n";
        }
    }

    edump.close();
    filenameStream.str("");

    // calculate condition numbers

    // calculate stability matrix
    // Ax=B --> x = linsolve(A,B)
    // t = linsolve(a,across)

    // double t[jnum*jtypes*jnum*jtypes];

    gsl_matrix_view a_gsl = gsl_matrix_view_array (a, jnum*jtypes, jnum*jtypes);

    for (i=0; i<(jnum*jtypes); i++)
    {
        for (j=0; j<(jnum*jtypes); j++)
        {
            across_i[j] = across[j*(jnum*jtypes)+i];
        }

        gsl_vector_view across_gsl_i = gsl_vector_view_array (across_i, (jnum*jtypes));

        gsl_vector *t_i = gsl_vector_alloc (jnum*jtypes);

        gsl_permutation * p = gsl_permutation_alloc (jnum*jtypes);

        gsl_linalg_LU_decomp (&a_gsl.matrix, p, &s);
        gsl_linalg_LU_solve (&a_gsl.matrix, p, &across_gsl_i.vector, t_i);

        filenameStream << "./t_stab.dat";
        filename = filenameStream.str();
        tdump.open(filename.c_str(), std::ios_base::app);

        for (j=0; j<(jnum*jtypes); j++)
        {
            tdump << gsl_vector_get (t_i, j) << "\n";
        }

        tdump.close();
        filenameStream.str("");

        for (j=0; j<(jnum*jtypes); j++)
        {
            t[j*(jnum*jtypes)+i] = gsl_vector_get (t_i, j);
        }

        gsl_permutation_free (p);
        gsl_vector_free (t_i);
    }

    // condition number:
    // int gsl_linalg_cholesky_rcond (const gsl_matrix * cholesky, double * rcond, gsl_vector * work)

    // eigensystem:

    gsl_matrix_view t_gsl = gsl_matrix_view_array (t, (jnum*jtypes), (jnum*jtypes));

    gsl_vector *s_svd = gsl_vector_alloc (jnum*jtypes);
    gsl_matrix *v_svd = gsl_matrix_alloc (jnum*jtypes, jnum*jtypes);
    gsl_vector *w_svd = gsl_vector_alloc (jnum*jtypes);

    //gsl_vector *eval = gsl_vector_alloc (jnum*jtypes);
    //gsl_matrix *evec = gsl_matrix_alloc (jnum*jtypes, jnum*jtypes);
    gsl_vector_complex *eval = gsl_vector_complex_alloc (jnum*jtypes);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (jnum*jtypes, jnum*jtypes);

    //gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (jnum*jtypes);
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (jnum*jtypes);

    //gsl_eigen_symmv (&t_gsl.matrix, eval, evec, w);
    gsl_eigen_nonsymmv (&t_gsl.matrix, eval, evec, w);

    gsl_linalg_SV_decomp (&t_gsl.matrix, v_svd, s_svd, w_svd);

    //gsl_eigen_symmv_free (w);
    gsl_eigen_nonsymmv_free (w);

    //gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
    gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    for (i=0; i<(jnum*jtypes); i++)
    {
        //double eval_i = gsl_vector_get (eval, i);
        //gsl_vector_view evec_i = gsl_matrix_column (evec, i);
        gsl_complex eval_i = gsl_vector_complex_get (eval, i);
        gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);

        filenameStream << "./eigvals.dat";
        filename = filenameStream.str();
        edump.open(filename.c_str(), std::ios_base::app);

        edump << GSL_REAL(eval_i) << "\t" << GSL_IMAG(eval_i) << "\n";

        edump.close();
        filenameStream.str("");

        filenameStream << "./exponent.dat";
        filename = filenameStream.str();
        edump.open(filename.c_str(), std::ios_base::app);

        //edump << log(GSL_REAL(eval_i))/log(b) << "\n";
        if (GSL_REAL(eval_i)>0.0)
        {
            edump << log(GSL_REAL(eval_i))/log(2) << "\n";
        }
        else
        {
            edump << -1000 << "\n";
        }

        edump.close();
        filenameStream.str("");

        filenameStream << "./eigvecs.dat";
        filename = filenameStream.str();
        edump.open(filename.c_str(), std::ios_base::app);

        for (j=0; j<(jnum*jtypes); j++)
        {
            gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
            edump << GSL_REAL(z) << "\t" << GSL_IMAG(z) << "\n";
        }

        edump.close();
        filenameStream.str("");

        filenameStream << "./singularvalues.dat";
        filename = filenameStream.str();
        edump.open(filename.c_str(), std::ios_base::app);

        for (j=0; j<(jnum*jtypes); j++)
        {
            edump << gsl_vector_get(s_svd, j) << "\n";
        }

        edump.close();
        filenameStream.str("");

        filenameStream << "./svd_v.dat";
        filename = filenameStream.str();
        edump.open(filename.c_str(), std::ios_base::app);

        for (j=0; j<(jnum*jtypes); j++)
        {
            for (k=0; k<(jnum*jtypes); k++)
            {
                // print out transpose(V)
                edump << gsl_matrix_get(v_svd, k, j) << "\n";
            }
        }

        edump.close();
        filenameStream.str("");

        gsl_vector_free (s_svd);
        gsl_matrix_free (v_svd);
        gsl_vector_free (w_svd);

        //gsl_vector_free (eval);
        //gsl_matrix_free (evec);
        gsl_vector_complex_free (eval);
        gsl_matrix_complex_free (evec);
    }

    filenameStream << "./jflow.dat";
    filename = filenameStream.str();
    edump.open(filename.c_str(), std::ios_base::app);

    for (i=0; i<(jtypes*jnum); i++)
    {
        edump << jvec[i] << "\n";
    }

    edump.close();
    filenameStream.str("");
}

critexp::~critexp()
{
}
