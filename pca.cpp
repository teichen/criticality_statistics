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
    int phase_pts;
    phase_pts = 1*6*11;

    int mu,c,lambda;

    int mu_vec[phase_pts];
    int c_vec[phase_pts];
    int lambda_vec[phase_pts];

    std::ostringstream filenameStream;
    std:string filename;

    ifstream fin;
    string line;

    std::ofstream tdump;
    std::ofstream edump;

    int i,j,k;
    int ii,jj,kk;

    int L[2];
    int npts[2];

    L[0] = 8;
    L[1] = 16;
    npts[0] = 20;
    npts[1] = 20;

    int nruns;
    //nruns = 10;
    nruns = 8; // memory constraints

    int runpt;
    double nave;

    jj = 0;
    for (mu=20; mu<=20; mu=mu+20)
    {
        for (c=50; c<=100; c=c+10)
        {
            for (lambda=50; lambda<=100; lambda=lambda+5)
            {
                mu_vec[jj] = mu;
                c_vec[jj] = c;
                lambda_vec[jj] = lambda;

                jj++;
            }
        }
    }

    double vt;

    // npts sampled at equilibrium

    double n[nruns*npts[1]*L[1]*L[1]*L[1]];
    double n_trans[nruns*npts[1]*L[1]*L[1]*L[1]];

    int s;

    for (jj=0; jj<phase_pts; jj++)
    {
        mu = mu_vec[jj];
        c = c_vec[jj];
        lambda = lambda_vec[jj];

        for (i=1; i<=nruns; i++)
        {
            for (runpt=0; runpt<npts[1]; runpt++)
            {
                // read in data

                filenameStream << "./L" << L[1] << "/" << mu << "/c" << c << "/" << lambda << "/" << i << "/netlib/n" << runpt << ".dat";
                filename = filenameStream.str();

                fin.open(filename.c_str());

                k = (i-1)*npts[1]*L[1]*L[1]*L[1] + runpt*L[1]*L[1]*L[1];
                while( getline(fin,line) )
                {
                    n[k] = atoi (line.c_str());
                    k++;
                }

                fin.close();
                filenameStream.str("");
            }
        }

        // shift relative to site means

        for (j=0; j<(L[1]*L[1]*L[1]); j++)
        {
            nave = 0.0;
            for (i=1; i<=nruns; i++)
            {
                for (runpt=0; runpt<npts[1]; runpt++)
                {
                    nave = nave + n[(i-1)*npts[1]*L[1]*L[1]*L[1]+runpt*L[1]*L[1]*L[1]+j];
                }
            }

            nave = nave/(double)(npts[1]*nruns);

            for (i=1; i<=nruns; i++)
            {
                for (runpt=0; runpt<npts[1]; runpt++)
                {
                    n[(i-1)*npts[1]*L[1]*L[1]*L[1]+runpt*L[1]*L[1]*L[1]+j] = n[(i-1)*npts[1]*L[1]*L[1]*L[1]+runpt*L[1]*L[1]*L[1]+j] - nave;
                }
            }
        }

        // calculate transpose of data

        for (j=0; j<(L[1]*L[1]*L[1]); j++)
        {
            for (i=1; i<=nruns; i++)
            {
                for (runpt=0; runpt<npts[1]; runpt++)
                {
                    n_trans[j*npts[1]*nruns+(i-1)*npts[1]+runpt] = n[(i-1)*npts[1]*L[1]*L[1]*L[1]+runpt*L[1]*L[1]*L[1]+j];
                }
            }
        }

        // PCA (SVD algorithm)

        // PCA of transpose of data

        gsl_matrix_view n_gsl = gsl_matrix_view_array (n_trans, L[1]*L[1]*L[1], npts[1]*nruns);

        gsl_vector *s_svd = gsl_vector_alloc (npts[1]*nruns);
        gsl_matrix *v_svd = gsl_matrix_alloc (npts[1]*nruns, npts[1]*nruns);
        gsl_vector *w_svd = gsl_vector_alloc (npts[1]*nruns);

        gsl_linalg_SV_decomp (&n_gsl.matrix, v_svd, s_svd, w_svd);

        filenameStream << "./pcvals.dat";
        filename = filenameStream.str();
        edump.open(filename.c_str(), std::ios_base::app);

        for (j=0; j<(npts[1]*nruns); j++)
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

        for (i=0; i<(npts[1]*nruns); i++)
        {
            for (j=0; j<L[1]; j++)
            {
                vt = 0.0;
                for (k=0; k<(npts[1]*nruns); k++)
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

}

pca::~pca()
{
}
