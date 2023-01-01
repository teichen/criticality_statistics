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
    int bind;

    int blevels,phase_pts,jtypes,jnum;
    blevels = 4;

/*//for((mu=20;mu<=40;mu=mu+20));
for((mu=20;mu<=20;mu=mu+20));
for((c=50;c<=100;c=c+10));
for((lambda=50;lambda<=100;lambda=lambda+5));
*/

    phase_pts = 1*6*11;
    jtypes = 5;
    jnum = 8;

    double jvec[blevels*phase_pts*jtypes*jnum];
    int jind;

    double nsum;
    int nswitch;

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
    nruns = 10;

    // npts sampled at equilibrium
    double na[nruns*npts[1]*jnum];
    double nb[nruns*npts[1]*jnum];
    double ua[nruns*npts[1]*jnum];
    double ub[nruns*npts[1]*jnum];
    double ma[nruns*npts[1]*jnum];
    double mb[nruns*npts[1]*jnum];
    double xa[nruns*npts[1]*jnum];
    double xb[nruns*npts[1]*jnum];
    double ca[nruns*npts[1]*jnum];
    double cb[nruns*npts[1]*jnum];

    // initialize jnum*2 matrices for mave,nave averages
    double nave[jnum*2];
    double uave[jnum*2];
    double mave[jnum*2];
    double xave[jnum*2];
    double cave[jnum*2];

    // initialize jnum*jnum matrices for mk,mj and nk,nj correlation functions
    double n_n[jnum*jnum];
    double n_u[jnum*jnum];
    double n_m[jnum*jnum];
    double n_x[jnum*jnum];
    double n_c[jnum*jnum];
    double u_n[jnum*jnum];
    double u_u[jnum*jnum];
    double u_m[jnum*jnum];
    double u_x[jnum*jnum];
    double u_c[jnum*jnum];
    double m_n[jnum*jnum];
    double m_u[jnum*jnum];
    double m_m[jnum*jnum];
    double m_x[jnum*jnum];
    double m_c[jnum*jnum];
    double x_n[jnum*jnum];
    double x_u[jnum*jnum];
    double x_m[jnum*jnum];
    double x_x[jnum*jnum];
    double x_c[jnum*jnum];
    double c_n[jnum*jnum];
    double c_u[jnum*jnum];
    double c_m[jnum*jnum];
    double c_x[jnum*jnum];
    double c_c[jnum*jnum];

    double a[jnum*jtypes*jnum*jtypes];
    double across[jnum*jtypes*jnum*jtypes];

    double t[jnum*jtypes*jnum*jtypes];

    double across_i[jnum*jtypes];
    int s;

    // read in data

    // collective distance from gel point:

    double d[phase_pts];

    filenameStream << "./d.dat";
    filename = filenameStream.str();

    fin.open(filename.c_str());

    j = 0;
    while( getline(fin,line) )
    {
        d[j] = atof (line.c_str());

        j++;
    }

    // minimum d[] value corresponds to gel point (lambda value)

    double j_gel[jtypes*jnum];

    int gp;
    i = 1; gp = 0;
    while (i<phase_pts)
    {
        if (d[i]<d[gp])
        {
            gp = i;
        }
        i++;
    }

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

    //lambda = 1.0*1.51*temp; // lattice gas scalar
    //mu = -3.0*lambda + (2.25e-4)*temp; // chemical potential
    //c = 60.0*pow(0.21,3)*temp;

    //lambda = 0.01*lambda_in*lambda;
    //c = 0.01*c_in*c;

    // K = 0.5, a = 4.2*temp/rhol

    mu = mu_vec[gp];
    c = c_vec[gp];
    lambda = lambda_vec[gp];

    double rhol;
    //rhol = 0.2734; // need the critical density here
    rhol = 0.005;

    j_gel[0] = (-3.0*(1.0*1.51*lambda*0.01)+(2.25e-4)); // mu/temp
    j_gel[jnum] = -(c*0.01*(60.0*pow(0.21,3))-2*0.5*(4.2/rhol)*pow(rhol,2)
                    +pow(0.21,3)*rhol/(double)(L[1]*L[1]*L[1])); // -(c-2*K*a*(rhol^2)+temp*(\ell^3/V)*rhol)/temp
    j_gel[2*jnum] = -(0.5*rhol); // -rhol/2
    j_gel[3*jnum] = 0.0;
    j_gel[4*jnum] = -3.0*(0.5*(4.2/rhol)*pow(rhol,2)/3.0); // -3.0*(K*a*pow(rhol,2)/3.0)/temp

    j_gel[1] = (1.0*1.51*lambda*0.01); // lambda/temp
    j_gel[jnum+1] = -(0.5*(4.2/rhol)*pow(rhol,2)/3.0); // -(K*a*pow(rhol,2)/3.0)/temp
    j_gel[2*jnum+1] = 0.0;
    j_gel[3*jnum+1] = 0.0;
    j_gel[4*jnum+1] = 0.0;

    for (i=2; i<jnum; i++)
    {
        j_gel[i] = 0.0;
        j_gel[jnum+i] = 0.0;
        j_gel[2*jnum+i] = 0.0;
        j_gel[3*jnum+i] = 0.0;
        j_gel[4*jnum+i] = 0.0;
    }

    for (i=0; i<(blevels*phase_pts*jtypes*jnum); i++)
    {
        jvec[i] = 0.0;
    }

    bind = 0;
    jind = 0;

    for (jj=0; jj<phase_pts; jj++)
    {
        mu = mu_vec[jj];
        c = c_vec[jj];
        lambda = lambda_vec[jj];

//lambda = 1.0*1.51*temp; // lattice gas scalar
//mu = -3.0*lambda + (2.25e-4)*temp; // chemical potential
//c = 60.0*pow(0.21,3)*temp;

//lambda = 0.01*lambda_in*lambda;
//c = 0.01*c_in*c;

// K = 0.5, a = 4.2*temp/rhol

        jvec[bind*phase_pts*jtypes*jnum+jind*jtypes*jnum] = (-3.0*(1.0*1.51*lambda*0.01)+(2.25e-4)); // mu/temp
        jvec[bind*phase_pts*jtypes*jnum+jind*jtypes*jnum+jnum] = -(c*0.01*(60.0*pow(0.21,3))-2*0.5*(4.2/rhol)*pow(rhol,2)
                                                                   +pow(0.21,3)*rhol/(double)(L[1]*L[1]*L[1])); // -(c-2*K*a*(rhol^2)+temp*(\ell^3/V)*rhol)/temp
        jvec[bind*phase_pts*jtypes*jnum+jind*jtypes*jnum+2*jnum] = -(0.5*rhol); // -rhol/2
        jvec[bind*phase_pts*jtypes*jnum+jind*jtypes*jnum+3*jnum] = 0.0;
        jvec[bind*phase_pts*jtypes*jnum+jind*jtypes*jnum+4*jnum] = -3.0*(0.5*(4.2/rhol)*pow(rhol,2)/3.0); // -3.0*(K*a*pow(rhol,2)/3.0)/temp

        jvec[bind*phase_pts*jtypes*jnum+jind*jtypes*jnum+1] = (1.0*1.51*lambda*0.01); // lambda/temp
        jvec[bind*phase_pts*jtypes*jnum+jind*jtypes*jnum+jnum+1] = -(0.5*(4.2/rhol)*pow(rhol,2)/3.0); // -(K*a*pow(rhol,2)/3.0)/temp
        jvec[bind*phase_pts*jtypes*jnum+jind*jtypes*jnum+2*jnum+1] = 0.0;
        jvec[bind*phase_pts*jtypes*jnum+jind*jtypes*jnum+3*jnum+1] = 0.0;
        jvec[bind*phase_pts*jtypes*jnum+jind*jtypes*jnum+4*jnum+1] = 0.0;

        jind++;
    }

    bind = 1;
    for (b=2; b<=8; b=b*2)
    {
        jind = 0;

        for (jj=0; jj<phase_pts; jj++)
        {
            mu = mu_vec[jj];
            c = c_vec[jj];
            lambda = lambda_vec[jj];

            for (i=1; i<=nruns; i++)
            {
                for (j=0; j<(npts[1]*jnum); j++)
                {
                    na[(i-1)*npts[1]*jnum+j] = 0.0;
                    ua[(i-1)*npts[1]*jnum+j] = 0.0;
                    ma[(i-1)*npts[1]*jnum+j] = 0.0;
                    xa[(i-1)*npts[1]*jnum+j] = 0.0;
                    ca[(i-1)*npts[1]*jnum+j] = 0.0;
                }
                for (j=0; j<(npts[1]*jnum); j++)
                {
                    nb[(i-1)*npts[1]*jnum+j] = 0.0;
                    ub[(i-1)*npts[1]*jnum+j] = 0.0;
                    mb[(i-1)*npts[1]*jnum+j] = 0.0;
                    xb[(i-1)*npts[1]*jnum+j] = 0.0;
                    cb[(i-1)*npts[1]*jnum+j] = 0.0;
                }

                filenameStream << "./L" << L[1] << "/" << mu << "/c" << c << "/" << lambda << "/" << i << "/nvec_b" << b/2 << ".dat";
                filename = filenameStream.str();

                fin.open(filename.c_str());

                j = (i-1)*npts[1]*jnum;
                while( getline(fin,line) )
                {
                    nb[j] = atof (line.c_str());

                    j++;
                }

                fin.close();
                filenameStream.str("");

                filenameStream << "./L" << L[1] << "/" << mu << "/c" << c << "/" << lambda << "/" << i << "/nvec_b" << b << ".dat";
                filename = filenameStream.str();

                fin.open(filename.c_str());

                j = (i-1)*npts[1]*jnum;
                while( getline(fin,line) )
                {
                    na[j] = atof (line.c_str());

                    j++;
                }

                fin.close();
                filenameStream.str("");

                filenameStream << "./L" << L[1] << "/" << mu << "/c" << c << "/" << lambda << "/" << i << "/uvec_b" << b/2 << ".dat";
                filename = filenameStream.str();

                fin.open(filename.c_str());

                j = (i-1)*npts[1]*jnum;
                while( getline(fin,line) )
                {
                    ub[j] = atof (line.c_str());
                    j++;
                }

                fin.close();
                filenameStream.str("");

                filenameStream << "./L" << L[1] << "/" << mu << "/c" << c << "/" << lambda << "/" << i << "/uvec_b" << b << ".dat";
                filename = filenameStream.str();

                fin.open(filename.c_str());

                j = (i-1)*npts[1]*jnum;
                while( getline(fin,line) )
                {
                    ua[j] = atof (line.c_str());
                    j++;
                }

                fin.close();
                filenameStream.str("");

                filenameStream << "./L" << L[1] << "/" << mu << "/c" << c << "/" << lambda << "/" << i << "/mvec_b" << b/2 << ".dat";
                filename = filenameStream.str();

                fin.open(filename.c_str());

                j = (i-1)*npts[1]*jnum;
                while( getline(fin,line) )
                {
                    mb[j] = atof (line.c_str());
                    j++;
                }

                fin.close();
                filenameStream.str("");

                filenameStream << "./L" << L[1] << "/" << mu << "/c" << c << "/" << lambda << "/" << i << "/mvec_b" << b << ".dat";
                filename = filenameStream.str();

                fin.open(filename.c_str());

                j = (i-1)*npts[1]*jnum;
                while( getline(fin,line) )
                {
                    ma[j] = atof (line.c_str());
                    j++;
                }

                fin.close();
                filenameStream.str("");

                filenameStream << "./L" << L[1] << "/" << mu << "/c" << c << "/" << lambda << "/" << i << "/xvec_b" << b/2 << ".dat";
                filename = filenameStream.str();

                fin.open(filename.c_str());

                j = (i-1)*npts[1]*jnum;
                while( getline(fin,line) )
                {
                    xb[j] = atof (line.c_str());
                    j++;
                }

                fin.close();
                filenameStream.str("");

                filenameStream << "./L" << L[1] << "/" << mu << "/c" << c << "/" << lambda << "/" << i << "/xvec_b" << b << ".dat";
                filename = filenameStream.str();

                fin.open(filename.c_str());

                j = (i-1)*npts[1]*jnum;
                while( getline(fin,line) )
                {
                    xa[j] = atof (line.c_str());
                    j++;
                }

                fin.close();
                filenameStream.str("");

                filenameStream << "./L" << L[1] << "/" << mu << "/c" << c << "/" << lambda << "/" << i << "/cvec_b" << b/2 << ".dat";
                filename = filenameStream.str();

                fin.open(filename.c_str());

                j = (i-1)*npts[1]*jnum;
                while( getline(fin,line) )
                {
                    cb[j] = atof (line.c_str());
                    j++;
                }

                fin.close();
                filenameStream.str("");

                filenameStream << "./L" << L[1] << "/" << mu << "/c" << c << "/" << lambda << "/" << i << "/cvec_b" << b << ".dat";
                filename = filenameStream.str();

                fin.open(filename.c_str());

                j = (i-1)*npts[1]*jnum;
                while( getline(fin,line) )
                {
                    ca[j] = atof (line.c_str());
                    j++;
                }

                fin.close();
                filenameStream.str("");

            }

            nswitch = 0;
            nsum = 0.0;

            for (j=0; j<(nruns*npts[1]*jnum); j++)
            {
                if (ua[j]>0.0)
                {
                    nsum = nsum+ua[j];
                }
            }

            if (nsum>0.0)
            {
                nswitch = 1;
            }

            // calculate averages
            for (i=0; i<jnum; i++)
            {
                nave[i*2+0] = 0.0;
                nave[i*2+1] = 0.0;
                uave[i*2+0] = 0.0;
                uave[i*2+1] = 0.0;
                mave[i*2+0] = 0.0;
                mave[i*2+1] = 0.0;
                xave[i*2+0] = 0.0;
                xave[i*2+1] = 0.0;
                cave[i*2+0] = 0.0;
                cave[i*2+1] = 0.0;

                for (j=0; j<(nruns*npts[1]); j++)
                {
                    nave[i*2+0] = nave[i*2+0] + na[j*jnum+i];
                    uave[i*2+0] = uave[i*2+0] + ua[j*jnum+i];
                    mave[i*2+0] = mave[i*2+0] + ma[j*jnum+i];
                    xave[i*2+0] = xave[i*2+0] + xa[j*jnum+i];
                    cave[i*2+0] = cave[i*2+0] + ca[j*jnum+i];
                }
                for (j=0; j<(nruns*npts[1]); j++)
                {
                    nave[i*2+1] = nave[i*2+1] + nb[j*jnum+i];
                    uave[i*2+1] = uave[i*2+1] + ub[j*jnum+i];
                    mave[i*2+1] = mave[i*2+1] + mb[j*jnum+i];
                    xave[i*2+1] = xave[i*2+1] + xb[j*jnum+i];
                    cave[i*2+1] = cave[i*2+1] + cb[j*jnum+i];
                }

                nave[i*2+0] = nave[i*2+0]/(double)(nruns*npts[1]);
                nave[i*2+1] = nave[i*2+1]/(double)(nruns*npts[1]);
                uave[i*2+0] = uave[i*2+0]/(double)(nruns*npts[1]);
                uave[i*2+1] = uave[i*2+1]/(double)(nruns*npts[1]);
                mave[i*2+0] = mave[i*2+0]/(double)(nruns*npts[1]);
                mave[i*2+1] = mave[i*2+1]/(double)(nruns*npts[1]);
                xave[i*2+0] = xave[i*2+0]/(double)(nruns*npts[1]);
                xave[i*2+1] = xave[i*2+1]/(double)(nruns*npts[1]);
                cave[i*2+0] = cave[i*2+0]/(double)(nruns*npts[1]);
                cave[i*2+1] = cave[i*2+1]/(double)(nruns*npts[1]);
            }

            // calculate correlation functions
            for (i=0; i<jnum; i++)
            {
                for (j=0; j<jnum; j++)
                {
                    n_n[i*jnum+j] = 0.0;
                    n_u[i*jnum+j] = 0.0;
                    n_m[i*jnum+j] = 0.0;
                    n_x[i*jnum+j] = 0.0;
                    n_c[i*jnum+j] = 0.0;
                    u_n[i*jnum+j] = 0.0;
                    u_u[i*jnum+j] = 0.0;
                    u_m[i*jnum+j] = 0.0;
                    u_x[i*jnum+j] = 0.0;
                    u_c[i*jnum+j] = 0.0;
                    m_n[i*jnum+j] = 0.0;
                    m_u[i*jnum+j] = 0.0;
                    m_m[i*jnum+j] = 0.0;
                    m_x[i*jnum+j] = 0.0;
                    m_c[i*jnum+j] = 0.0;
                    x_n[i*jnum+j] = 0.0;
                    x_u[i*jnum+j] = 0.0;
                    x_m[i*jnum+j] = 0.0;
                    x_x[i*jnum+j] = 0.0;
                    x_c[i*jnum+j] = 0.0;
                    c_n[i*jnum+j] = 0.0;
                    c_u[i*jnum+j] = 0.0;
                    c_m[i*jnum+j] = 0.0;
                    c_x[i*jnum+j] = 0.0;
                    c_c[i*jnum+j] = 0.0;

                    for (k=0; k<(nruns*npts[1]); k++)
                    {
                        n_n[i*jnum+j] = n_n[i*jnum+j] + na[k*jnum+i]*na[k*jnum+j];
                        n_u[i*jnum+j] = n_u[i*jnum+j] + na[k*jnum+i]*ua[k*jnum+j];
                        n_m[i*jnum+j] = n_m[i*jnum+j] + na[k*jnum+i]*ma[k*jnum+j];
                        n_x[i*jnum+j] = n_x[i*jnum+j] + na[k*jnum+i]*xa[k*jnum+j];
                        n_c[i*jnum+j] = n_c[i*jnum+j] + na[k*jnum+i]*ca[k*jnum+j];
                        u_n[i*jnum+j] = u_n[i*jnum+j] + ua[k*jnum+i]*na[k*jnum+j];
                        u_u[i*jnum+j] = u_u[i*jnum+j] + ua[k*jnum+i]*ua[k*jnum+j];
                        u_m[i*jnum+j] = u_m[i*jnum+j] + ua[k*jnum+i]*ma[k*jnum+j];
                        u_x[i*jnum+j] = u_x[i*jnum+j] + ua[k*jnum+i]*xa[k*jnum+j];
                        u_c[i*jnum+j] = u_c[i*jnum+j] + ua[k*jnum+i]*ca[k*jnum+j];
                        m_n[i*jnum+j] = m_n[i*jnum+j] + ma[k*jnum+i]*na[k*jnum+j];
                        m_u[i*jnum+j] = m_u[i*jnum+j] + ma[k*jnum+i]*ua[k*jnum+j];
                        m_m[i*jnum+j] = m_m[i*jnum+j] + ma[k*jnum+i]*ma[k*jnum+j];
                        m_x[i*jnum+j] = m_x[i*jnum+j] + ma[k*jnum+i]*xa[k*jnum+j];
                        m_c[i*jnum+j] = m_c[i*jnum+j] + ma[k*jnum+i]*ca[k*jnum+j];
                        x_n[i*jnum+j] = x_n[i*jnum+j] + xa[k*jnum+i]*na[k*jnum+j];
                        x_u[i*jnum+j] = x_u[i*jnum+j] + xa[k*jnum+i]*ua[k*jnum+j];
                        x_m[i*jnum+j] = x_m[i*jnum+j] + xa[k*jnum+i]*ma[k*jnum+j];
                        x_x[i*jnum+j] = x_x[i*jnum+j] + xa[k*jnum+i]*xa[k*jnum+j];
                        x_c[i*jnum+j] = x_c[i*jnum+j] + xa[k*jnum+i]*ca[k*jnum+j];
                        c_n[i*jnum+j] = c_n[i*jnum+j] + ca[k*jnum+i]*na[k*jnum+j];
                        c_u[i*jnum+j] = c_u[i*jnum+j] + ca[k*jnum+i]*ua[k*jnum+j];
                        c_m[i*jnum+j] = c_m[i*jnum+j] + ca[k*jnum+i]*ma[k*jnum+j];
                        c_x[i*jnum+j] = c_x[i*jnum+j] + ca[k*jnum+i]*xa[k*jnum+j];
                        c_c[i*jnum+j] = c_c[i*jnum+j] + ca[k*jnum+i]*ca[k*jnum+j];
                    }
                    n_n[i*jnum+j] = n_n[i*jnum+j]/(double)(nruns*npts[1]);
                    n_u[i*jnum+j] = n_u[i*jnum+j]/(double)(nruns*npts[1]);
                    n_m[i*jnum+j] = n_m[i*jnum+j]/(double)(nruns*npts[1]);
                    n_x[i*jnum+j] = n_x[i*jnum+j]/(double)(nruns*npts[1]);
                    n_c[i*jnum+j] = n_c[i*jnum+j]/(double)(nruns*npts[1]);
                    u_n[i*jnum+j] = u_n[i*jnum+j]/(double)(nruns*npts[1]);
                    u_u[i*jnum+j] = u_u[i*jnum+j]/(double)(nruns*npts[1]);
                    u_m[i*jnum+j] = u_m[i*jnum+j]/(double)(nruns*npts[1]);
                    u_x[i*jnum+j] = u_x[i*jnum+j]/(double)(nruns*npts[1]);
                    u_c[i*jnum+j] = u_c[i*jnum+j]/(double)(nruns*npts[1]);
                    m_n[i*jnum+j] = m_n[i*jnum+j]/(double)(nruns*npts[1]);
                    m_u[i*jnum+j] = m_u[i*jnum+j]/(double)(nruns*npts[1]);
                    m_m[i*jnum+j] = m_m[i*jnum+j]/(double)(nruns*npts[1]);
                    m_x[i*jnum+j] = m_x[i*jnum+j]/(double)(nruns*npts[1]);
                    m_c[i*jnum+j] = m_c[i*jnum+j]/(double)(nruns*npts[1]);
                    x_n[i*jnum+j] = x_n[i*jnum+j]/(double)(nruns*npts[1]);
                    x_u[i*jnum+j] = x_u[i*jnum+j]/(double)(nruns*npts[1]);
                    x_m[i*jnum+j] = x_m[i*jnum+j]/(double)(nruns*npts[1]);
                    x_x[i*jnum+j] = x_x[i*jnum+j]/(double)(nruns*npts[1]);
                    x_c[i*jnum+j] = x_c[i*jnum+j]/(double)(nruns*npts[1]);
                    c_n[i*jnum+j] = c_n[i*jnum+j]/(double)(nruns*npts[1]);
                    c_u[i*jnum+j] = c_u[i*jnum+j]/(double)(nruns*npts[1]);
                    c_m[i*jnum+j] = c_m[i*jnum+j]/(double)(nruns*npts[1]);
                    c_x[i*jnum+j] = c_x[i*jnum+j]/(double)(nruns*npts[1]);
                    c_c[i*jnum+j] = c_c[i*jnum+j]/(double)(nruns*npts[1]);

                    n_n[i*jnum+j] = n_n[i*jnum+j] - nave[i*2+1]*nave[j*2+1];
                    n_u[i*jnum+j] = n_u[i*jnum+j] - nave[i*2+1]*uave[j*2+1];
                    n_m[i*jnum+j] = n_m[i*jnum+j] - nave[i*2+1]*mave[j*2+1];
                    n_x[i*jnum+j] = n_x[i*jnum+j] - nave[i*2+1]*xave[j*2+1];
                    n_c[i*jnum+j] = n_c[i*jnum+j] - nave[i*2+1]*cave[j*2+1];
                    u_n[i*jnum+j] = u_n[i*jnum+j] - uave[i*2+1]*nave[j*2+1];
                    u_u[i*jnum+j] = u_u[i*jnum+j] - uave[i*2+1]*uave[j*2+1];
                    u_m[i*jnum+j] = u_m[i*jnum+j] - uave[i*2+1]*mave[j*2+1];
                    u_x[i*jnum+j] = u_x[i*jnum+j] - uave[i*2+1]*xave[j*2+1];
                    u_c[i*jnum+j] = u_c[i*jnum+j] - uave[i*2+1]*cave[j*2+1];
                    m_n[i*jnum+j] = m_n[i*jnum+j] - mave[i*2+1]*nave[j*2+1];
                    m_u[i*jnum+j] = m_u[i*jnum+j] - mave[i*2+1]*uave[j*2+1];
                    m_m[i*jnum+j] = m_m[i*jnum+j] - mave[i*2+1]*mave[j*2+1];
                    m_x[i*jnum+j] = m_x[i*jnum+j] - mave[i*2+1]*xave[j*2+1];
                    m_c[i*jnum+j] = m_c[i*jnum+j] - mave[i*2+1]*cave[j*2+1];
                    x_n[i*jnum+j] = x_n[i*jnum+j] - xave[i*2+1]*nave[j*2+1];
                    x_u[i*jnum+j] = x_u[i*jnum+j] - xave[i*2+1]*uave[j*2+1];
                    x_m[i*jnum+j] = x_m[i*jnum+j] - xave[i*2+1]*mave[j*2+1];
                    x_x[i*jnum+j] = x_x[i*jnum+j] - xave[i*2+1]*xave[j*2+1];
                    x_c[i*jnum+j] = x_c[i*jnum+j] - xave[i*2+1]*cave[j*2+1];
                    c_n[i*jnum+j] = c_n[i*jnum+j] - cave[i*2+1]*nave[j*2+1];
                    c_u[i*jnum+j] = c_u[i*jnum+j] - cave[i*2+1]*uave[j*2+1];
                    c_m[i*jnum+j] = c_m[i*jnum+j] - cave[i*2+1]*mave[j*2+1];
                    c_x[i*jnum+j] = c_x[i*jnum+j] - cave[i*2+1]*xave[j*2+1];
                    c_c[i*jnum+j] = c_c[i*jnum+j] - cave[i*2+1]*cave[j*2+1];
                }
            }

            for (i=0; i<jnum; i++)
            {
                for (j=0; j<jnum; j++)
                {
                    a[i*jnum*jtypes+j] = n_n[i*jnum+j];
                    a[i*jnum*jtypes+j+jnum] = n_u[i*jnum+j];
                    a[i*jnum*jtypes+j+jnum*2] = n_m[i*jnum+j];
                    a[i*jnum*jtypes+j+jnum*3] = n_x[i*jnum+j];
                    a[i*jnum*jtypes+j+jnum*4] = n_c[i*jnum+j];
                    a[(i+jnum)*jnum*jtypes+j] = u_n[i*jnum+j];
                    a[(i+jnum)*jnum*jtypes+j+jnum] = u_u[i*jnum+j];
                    a[(i+jnum)*jnum*jtypes+j+jnum*2] = u_m[i*jnum+j];
                    a[(i+jnum)*jnum*jtypes+j+jnum*3] = u_x[i*jnum+j];
                    a[(i+jnum)*jnum*jtypes+j+jnum*4] = u_c[i*jnum+j];
                    a[(i+jnum*2)*jnum*jtypes+j] = m_n[i*jnum+j];
                    a[(i+jnum*2)*jnum*jtypes+j+jnum] = m_u[i*jnum+j];
                    a[(i+jnum*2)*jnum*jtypes+j+jnum*2] = m_m[i*jnum+j];
                    a[(i+jnum*2)*jnum*jtypes+j+jnum*3] = m_x[i*jnum+j];
                    a[(i+jnum*2)*jnum*jtypes+j+jnum*4] = m_c[i*jnum+j];
                    a[(i+jnum*3)*jnum*jtypes+j] = x_n[i*jnum+j];
                    a[(i+jnum*3)*jnum*jtypes+j+jnum] = x_u[i*jnum+j];
                    a[(i+jnum*3)*jnum*jtypes+j+jnum*2] = x_m[i*jnum+j];
                    a[(i+jnum*3)*jnum*jtypes+j+jnum*3] = x_x[i*jnum+j];
                    a[(i+jnum*3)*jnum*jtypes+j+jnum*4] = x_c[i*jnum+j];
                    a[(i+jnum*4)*jnum*jtypes+j] = c_n[i*jnum+j];
                    a[(i+jnum*4)*jnum*jtypes+j+jnum] = c_u[i*jnum+j];
                    a[(i+jnum*4)*jnum*jtypes+j+jnum*2] = c_m[i*jnum+j];
                    a[(i+jnum*4)*jnum*jtypes+j+jnum*3] = c_x[i*jnum+j];
                    a[(i+jnum*4)*jnum*jtypes+j+jnum*4] = c_c[i*jnum+j];
                }
            }

            for (i=0; i<jnum; i++)
            {
                for (j=0; j<jnum; j++)
                {
                    n_n[i*jnum+j] = 0.0;
                    n_u[i*jnum+j] = 0.0;
                    n_m[i*jnum+j] = 0.0;
                    n_x[i*jnum+j] = 0.0;
                    n_c[i*jnum+j] = 0.0;
                    u_n[i*jnum+j] = 0.0;
                    u_u[i*jnum+j] = 0.0;
                    u_m[i*jnum+j] = 0.0;
                    u_x[i*jnum+j] = 0.0;
                    u_c[i*jnum+j] = 0.0;
                    m_n[i*jnum+j] = 0.0;
                    m_u[i*jnum+j] = 0.0;
                    m_m[i*jnum+j] = 0.0;
                    m_x[i*jnum+j] = 0.0;
                    m_c[i*jnum+j] = 0.0;
                    x_n[i*jnum+j] = 0.0;
                    x_u[i*jnum+j] = 0.0;
                    x_m[i*jnum+j] = 0.0;
                    x_x[i*jnum+j] = 0.0;
                    x_c[i*jnum+j] = 0.0;
                    c_n[i*jnum+j] = 0.0;
                    c_u[i*jnum+j] = 0.0;
                    c_m[i*jnum+j] = 0.0;
                    c_x[i*jnum+j] = 0.0;
                    c_c[i*jnum+j] = 0.0;

                    for (k=0; k<(nruns*npts[1]); k++)
                    {
                        n_n[i*jnum+j] = n_n[i*jnum+j] + na[k*jnum+i]*nb[k*jnum+j];
                        n_u[i*jnum+j] = n_u[i*jnum+j] + na[k*jnum+i]*ub[k*jnum+j];
                        n_m[i*jnum+j] = n_m[i*jnum+j] + na[k*jnum+i]*mb[k*jnum+j];
                        n_x[i*jnum+j] = n_x[i*jnum+j] + na[k*jnum+i]*xb[k*jnum+j];
                        n_c[i*jnum+j] = n_c[i*jnum+j] + na[k*jnum+i]*cb[k*jnum+j];
                        u_n[i*jnum+j] = u_n[i*jnum+j] + ua[k*jnum+i]*nb[k*jnum+j];
                        u_u[i*jnum+j] = u_u[i*jnum+j] + ua[k*jnum+i]*ub[k*jnum+j];
                        u_m[i*jnum+j] = u_m[i*jnum+j] + ua[k*jnum+i]*mb[k*jnum+j];
                        u_x[i*jnum+j] = u_x[i*jnum+j] + ua[k*jnum+i]*xb[k*jnum+j];
                        u_c[i*jnum+j] = u_c[i*jnum+j] + ua[k*jnum+i]*cb[k*jnum+j];
                        m_n[i*jnum+j] = m_n[i*jnum+j] + ma[k*jnum+i]*nb[k*jnum+j];
                        m_u[i*jnum+j] = m_u[i*jnum+j] + ma[k*jnum+i]*ub[k*jnum+j];
                        m_m[i*jnum+j] = m_m[i*jnum+j] + ma[k*jnum+i]*mb[k*jnum+j];
                        m_x[i*jnum+j] = m_x[i*jnum+j] + ma[k*jnum+i]*xb[k*jnum+j];
                        m_c[i*jnum+j] = m_c[i*jnum+j] + ma[k*jnum+i]*cb[k*jnum+j];
                        x_n[i*jnum+j] = x_n[i*jnum+j] + xa[k*jnum+i]*nb[k*jnum+j];
                        x_u[i*jnum+j] = x_u[i*jnum+j] + xa[k*jnum+i]*ub[k*jnum+j];
                        x_m[i*jnum+j] = x_m[i*jnum+j] + xa[k*jnum+i]*mb[k*jnum+j];
                        x_x[i*jnum+j] = x_x[i*jnum+j] + xa[k*jnum+i]*xb[k*jnum+j];
                        x_c[i*jnum+j] = x_c[i*jnum+j] + xa[k*jnum+i]*cb[k*jnum+j];
                        c_n[i*jnum+j] = c_n[i*jnum+j] + ca[k*jnum+i]*nb[k*jnum+j];
                        c_u[i*jnum+j] = c_u[i*jnum+j] + ca[k*jnum+i]*ub[k*jnum+j];
                        c_m[i*jnum+j] = c_m[i*jnum+j] + ca[k*jnum+i]*mb[k*jnum+j];
                        c_x[i*jnum+j] = c_x[i*jnum+j] + ca[k*jnum+i]*xb[k*jnum+j];
                        c_c[i*jnum+j] = c_c[i*jnum+j] + ca[k*jnum+i]*cb[k*jnum+j];
                    }
                    n_n[i*jnum+j] = n_n[i*jnum+j]/(double)(nruns*npts[1]);
                    n_u[i*jnum+j] = n_u[i*jnum+j]/(double)(nruns*npts[1]);
                    n_m[i*jnum+j] = n_m[i*jnum+j]/(double)(nruns*npts[1]);
                    n_x[i*jnum+j] = n_x[i*jnum+j]/(double)(nruns*npts[1]);
                    n_c[i*jnum+j] = n_c[i*jnum+j]/(double)(nruns*npts[1]);
                    u_n[i*jnum+j] = u_n[i*jnum+j]/(double)(nruns*npts[1]);
                    u_u[i*jnum+j] = u_u[i*jnum+j]/(double)(nruns*npts[1]);
                    u_m[i*jnum+j] = u_m[i*jnum+j]/(double)(nruns*npts[1]);
                    u_x[i*jnum+j] = u_x[i*jnum+j]/(double)(nruns*npts[1]);
                    u_c[i*jnum+j] = u_c[i*jnum+j]/(double)(nruns*npts[1]);
                    m_n[i*jnum+j] = m_n[i*jnum+j]/(double)(nruns*npts[1]);
                    m_u[i*jnum+j] = m_u[i*jnum+j]/(double)(nruns*npts[1]);
                    m_m[i*jnum+j] = m_m[i*jnum+j]/(double)(nruns*npts[1]);
                    m_x[i*jnum+j] = m_x[i*jnum+j]/(double)(nruns*npts[1]);
                    m_c[i*jnum+j] = m_c[i*jnum+j]/(double)(nruns*npts[1]);
                    x_n[i*jnum+j] = x_n[i*jnum+j]/(double)(nruns*npts[1]);
                    x_u[i*jnum+j] = x_u[i*jnum+j]/(double)(nruns*npts[1]);
                    x_m[i*jnum+j] = x_m[i*jnum+j]/(double)(nruns*npts[1]);
                    x_x[i*jnum+j] = x_x[i*jnum+j]/(double)(nruns*npts[1]);
                    x_c[i*jnum+j] = x_c[i*jnum+j]/(double)(nruns*npts[1]);
                    c_n[i*jnum+j] = c_n[i*jnum+j]/(double)(nruns*npts[1]);
                    c_u[i*jnum+j] = c_u[i*jnum+j]/(double)(nruns*npts[1]);
                    c_m[i*jnum+j] = c_m[i*jnum+j]/(double)(nruns*npts[1]);
                    c_x[i*jnum+j] = c_x[i*jnum+j]/(double)(nruns*npts[1]);
                    c_c[i*jnum+j] = c_c[i*jnum+j]/(double)(nruns*npts[1]);

                    n_n[i*jnum+j] = n_n[i*jnum+j] - nave[i*2+1]*nave[j*2+0];
                    n_u[i*jnum+j] = n_u[i*jnum+j] - nave[i*2+1]*uave[j*2+0];
                    n_m[i*jnum+j] = n_m[i*jnum+j] - nave[i*2+1]*mave[j*2+0];
                    n_x[i*jnum+j] = n_x[i*jnum+j] - nave[i*2+1]*xave[j*2+0];
                    n_c[i*jnum+j] = n_c[i*jnum+j] - nave[i*2+1]*cave[j*2+0];
                    u_n[i*jnum+j] = u_n[i*jnum+j] - uave[i*2+1]*nave[j*2+0];
                    u_u[i*jnum+j] = u_u[i*jnum+j] - uave[i*2+1]*uave[j*2+0];
                    u_m[i*jnum+j] = u_m[i*jnum+j] - uave[i*2+1]*mave[j*2+0];
                    u_x[i*jnum+j] = u_x[i*jnum+j] - uave[i*2+1]*xave[j*2+0];
                    u_c[i*jnum+j] = u_c[i*jnum+j] - uave[i*2+1]*cave[j*2+0];
                    m_n[i*jnum+j] = m_n[i*jnum+j] - mave[i*2+1]*nave[j*2+0];
                    m_u[i*jnum+j] = m_u[i*jnum+j] - mave[i*2+1]*uave[j*2+0];
                    m_m[i*jnum+j] = m_m[i*jnum+j] - mave[i*2+1]*mave[j*2+0];
                    m_x[i*jnum+j] = m_x[i*jnum+j] - mave[i*2+1]*xave[j*2+0];
                    m_c[i*jnum+j] = m_c[i*jnum+j] - mave[i*2+1]*cave[j*2+0];
                    x_n[i*jnum+j] = x_n[i*jnum+j] - xave[i*2+1]*nave[j*2+0];
                    x_u[i*jnum+j] = x_u[i*jnum+j] - xave[i*2+1]*uave[j*2+0];
                    x_m[i*jnum+j] = x_m[i*jnum+j] - xave[i*2+1]*mave[j*2+0];
                    x_x[i*jnum+j] = x_x[i*jnum+j] - xave[i*2+1]*xave[j*2+0];
                    x_c[i*jnum+j] = x_c[i*jnum+j] - xave[i*2+1]*cave[j*2+0];
                    c_n[i*jnum+j] = c_n[i*jnum+j] - cave[i*2+1]*nave[j*2+0];
                    c_u[i*jnum+j] = c_u[i*jnum+j] - cave[i*2+1]*uave[j*2+0];
                    c_m[i*jnum+j] = c_m[i*jnum+j] - cave[i*2+1]*mave[j*2+0];
                    c_x[i*jnum+j] = c_x[i*jnum+j] - cave[i*2+1]*xave[j*2+0];
                    c_c[i*jnum+j] = c_c[i*jnum+j] - cave[i*2+1]*cave[j*2+0];
                }
            }

            for (i=0; i<jnum; i++)
            {
                for (j=0; j<jnum; j++)
                {
                    across[i*jnum*jtypes+j] = n_n[i*jnum+j];
                    across[i*jnum*jtypes+j+jnum] = n_u[i*jnum+j];
                    across[i*jnum*jtypes+j+jnum*2] = n_m[i*jnum+j];
                    across[i*jnum*jtypes+j+jnum*3] = n_x[i*jnum+j];
                    across[i*jnum*jtypes+j+jnum*4] = n_c[i*jnum+j];
                    across[(i+jnum)*jnum*jtypes+j] = u_n[i*jnum+j];
                    across[(i+jnum)*jnum*jtypes+j+jnum] = u_u[i*jnum+j];
                    across[(i+jnum)*jnum*jtypes+j+jnum*2] = u_m[i*jnum+j];
                    across[(i+jnum)*jnum*jtypes+j+jnum*3] = u_x[i*jnum+j];
                    across[(i+jnum)*jnum*jtypes+j+jnum*4] = u_c[i*jnum+j];
                    across[(i+jnum*2)*jnum*jtypes+j] = m_n[i*jnum+j];
                    across[(i+jnum*2)*jnum*jtypes+j+jnum] = m_u[i*jnum+j];
                    across[(i+jnum*2)*jnum*jtypes+j+jnum*2] = m_m[i*jnum+j];
                    across[(i+jnum*2)*jnum*jtypes+j+jnum*3] = m_x[i*jnum+j];
                    across[(i+jnum*2)*jnum*jtypes+j+jnum*4] = m_c[i*jnum+j];
                    across[(i+jnum*3)*jnum*jtypes+j] = x_n[i*jnum+j];
                    across[(i+jnum*3)*jnum*jtypes+j+jnum] = x_u[i*jnum+j];
                    across[(i+jnum*3)*jnum*jtypes+j+jnum*2] = x_m[i*jnum+j];
                    across[(i+jnum*3)*jnum*jtypes+j+jnum*3] = x_x[i*jnum+j];
                    across[(i+jnum*3)*jnum*jtypes+j+jnum*4] = x_c[i*jnum+j];
                    across[(i+jnum*4)*jnum*jtypes+j] = c_n[i*jnum+j];
                    across[(i+jnum*4)*jnum*jtypes+j+jnum] = c_u[i*jnum+j];
                    across[(i+jnum*4)*jnum*jtypes+j+jnum*2] = c_m[i*jnum+j];
                    across[(i+jnum*4)*jnum*jtypes+j+jnum*3] = c_x[i*jnum+j];
                    across[(i+jnum*4)*jnum*jtypes+j+jnum*4] = c_c[i*jnum+j];
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

            // temporarily ignore c-field

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

                if (nswitch==1)
                {
                    gsl_linalg_LU_decomp (&a_gsl.matrix, p, &s);
                    gsl_linalg_LU_solve (&a_gsl.matrix, p, &across_gsl_i.vector, t_i);
                }

                filenameStream << "./t_stab.dat";
                filename = filenameStream.str();
                tdump.open(filename.c_str(), std::ios_base::app);

                if (nswitch==1)
                {
                    for (j=0; j<(jnum*jtypes); j++)
                    {
                        tdump << gsl_vector_get (t_i, j) << "\n";
                    }
                }
                else
                {
                    for (j=0; j<(jnum*jtypes); j++)
                    {
                        tdump << 0.0 << "\n";
                    }
                }

                tdump.close();
                filenameStream.str("");

                if (nswitch==1)
                {
                    for (j=0; j<(jnum*jtypes); j++)
                    {
                        t[j*(jnum*jtypes)+i] = gsl_vector_get (t_i, j);
                    }
                }
                else
                {
                    for (j=0; j<(jnum*jtypes); j++)
                    {
                        t[j*(jnum*jtypes)+i] = 0.0;
                    }
                }

                gsl_permutation_free (p);
                gsl_vector_free (t_i);
            }

            for (i=0; i<(jnum*jtypes); i++)
            {
                jvec[bind*phase_pts*jtypes*jnum+jind*jtypes*jnum+i] = j_gel[i];
                for (j=0; j<(jnum*jtypes); j++)
                {
                    jvec[bind*phase_pts*jtypes*jnum+jind*jtypes*jnum+i] = jvec[bind*phase_pts*jtypes*jnum+jind*jtypes*jnum+i] + t[i*(jnum*jtypes)+j]*(jvec[(bind-1)*phase_pts*jtypes*jnum+jind*jtypes*jnum+j]-j_gel[j]);
                }
            }

            // condition number:
            // int gsl_linalg_cholesky_rcond (const gsl_matrix * cholesky, double * rcond, gsl_vector * work)

            // eigensystem:

            if (nswitch==1)
            {
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

                }

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
            else
            {
                for (i=0; i<(jnum*jtypes); i++)
                {
                    filenameStream << "./eigvals.dat";
                    filename = filenameStream.str();
                    edump.open(filename.c_str(), std::ios_base::app);

                    edump << 0.0 << "\t" << 0.0 << "\n";
                    edump.close();
                    filenameStream.str("");

                    filenameStream << "./exponent.dat";
                    filename = filenameStream.str();
                    edump.open(filename.c_str(), std::ios_base::app);

                    edump << -1000 << "\n";
                    edump.close();
                    filenameStream.str("");

                    filenameStream << "./eigvecs.dat";
                    filename = filenameStream.str();
                    edump.open(filename.c_str(), std::ios_base::app);

                    for (j=0; j<(jnum*jtypes); j++)
                    {
                        edump << 0.0 << "\t" << 0.0 << "\n";
                    }

                    edump.close();
                    filenameStream.str("");
                }

                filenameStream << "./singularvalues.dat";
                filename = filenameStream.str();
                edump.open(filename.c_str(), std::ios_base::app);

                for (j=0; j<(jnum*jtypes); j++)
                {
                    edump << -1000 << "\n";
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
                        edump << 0.0 << "\n";
                    }
                }

                edump.close();
                filenameStream.str("");

            }

            jind++;
        }
        bind++;
    }

    filenameStream << "./jflow.dat";
    filename = filenameStream.str();
    edump.open(filename.c_str(), std::ios_base::app);

    for (i=0; i<(blevels*phase_pts*jtypes*jnum); i++)
    {
        edump << jvec[i] << "\n";
    }

    edump.close();
    filenameStream.str("");


}

critexp::~critexp()
{
}
