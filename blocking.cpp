/*
    The blocking class calculates blocking statistics of a 
    binary discretized field.
*/
#include "blocking.h"
#include <cstdlib>
#include <stdlib.h>
#include <fstream>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <sstream>
using std::ifstream;
#include <cstring>

#include <sys/time.h>
#include <unistd.h>

#include <time.h>
#include <math.h>
#include <cmath>

#include <omp.h>

const double PI = 3.1415926535897932384626433832795028841971693;

using namespace std;

blocking::blocking()
{
}

blocking::blocking(int b)
{
    L = lattice.L; // length of lattice (number of sites)
    Lfine = lattice.Lfine; // fine grid

    dim = lattice.dim; // dimensionality of lattice

    nL = pow(L,dim);
    nLfine = pow(Lfine,dim);

    int runpts;
    runpts = 20;

    nvols = 100;

    srand (time(NULL)+b); // re-initialize random seed

    nL = pow(L,dim);
    nLfine = pow(Lfine,dim);

    // Fluctuation method of calculating elasticity

    std::ostringstream filenameStream;
    std:string filename;

    ifstream fin;
    string line;

    std::ofstream ndump;

    int i,j,k;
    int ib,jb,kb;

    // read in data
    // runtype tells you which block of trajectory data to calculate on

    int n[nL];

    double* via;
    via = (double*) calloc (nL*nLfine, sizeof(double));

    for (i=0; i<(int)(nL*nLfine); i++)
    {
        via[i] = 0.0;
    }

    double v[nL];
    double c[nL];

    // coarse-grained n and v (averages per cell)
    double ni[(int)pow(L/b,dim)];
    double vi[(int)pow(L/b,dim)];
    double ci[(int)pow(L/b,dim)];

    // collective variables for each MC config
    double n_1; // sum\, n_i
    double n0; // sum_NN\, n_i * n_i+1 (i=x,y,z)
    double n1; // sum_{diagonal in plane}\, n_i * n_j
    // e.g. n_(0,0,0) * n_(1,1,0)
    double n2; // sum_{cubic diagonal}\, n_i * n_j
    // e.g. n_(0,0,0) * n_(1,1,1)
    double n3; // sum_{principal planes}\, n_i * n_j * n_k * n_l
    // e.g. n_(0,0,0) * n_(0,1,0) * n_(1,0,0) * n_(1,1,0)
    double n4; // sum_{diagonal planes}\, n_i * n_j * n_k * n_l
    // e.g. n_(0,0,0) * n_(1,0,0) * n_(0,1,1) * n_(1,1,1)
    double n5; // sum_{tetrahedral vertices}\, n_i * n_j * n_k * n_l
    // e.g. n_(0,0,0) * n_(1,0,1) * n_(0,1,1) * n_(1,1,0)
    double n6; // sum_NNN\, n_i * n_i+2 (i=x,y,z)

    double u_1; // sum\, n_i*v_i
    double u0; // sum_NN\, n_i*v_i * n_i+1 (i=x,y,z)
    double u1; // sum_{diagonal in plane}\, n_i*v_i * n_j
    // e.g. n_(0,0,0) * n_(1,1,0)
    double u2; // sum_{cubic diagonal}\, n_i*v_i * n_j
    // e.g. n_(0,0,0) * n_(1,1,1)
    double u3; // sum_{principal planes}\, n_i*v_i * n_j * n_k*v_k * n_l
    // e.g. n_(0,0,0) * n_(0,1,0) * n_(1,0,0) * n_(1,1,0)
    double u4; // sum_{diagonal planes}\, n_i*v_i * n_j * n_k*v_k * n_l
    // e.g. n_(0,0,0) * n_(1,0,0) * n_(0,1,1) * n_(1,1,1)
    double u5; // sum_{tetrahedral vertices}\, n_i*v_i * n_j * n_k*v_k * n_l
    // e.g. n_(0,0,0) * n_(1,0,1) * n_(0,1,1) * n_(1,1,0)
    double u6; // sum_NNN\, n_i*v_i * n_i+2 (i=x,y,z)

    double m_1; // sum\, n_i^2
    double m0; // sum_NN\, n_i^2 * n_i+1 (i=x,y,z)
    double m1; // sum_{diagonal in plane}\, n_i^2 * n_j
    // e.g. n_(0,0,0) * n_(1,1,0)
    double m2; // sum_{cubic diagonal}\, n_i^2 * n_j
    // e.g. n_(0,0,0) * n_(1,1,1)
    double m3; // sum_{principal planes}\, n_i^2 * n_j * n_k^2 * n_l
    // e.g. n_(0,0,0) * n_(0,1,0) * n_(1,0,0) * n_(1,1,0)
    double m4; // sum_{diagonal planes}\, n_i^2 * n_j * n_k^2 * n_l
    // e.g. n_(0,0,0) * n_(1,0,0) * n_(0,1,1) * n_(1,1,1)
    double m5; // sum_{tetrahedral vertices}\, n_i^2 * n_j * n_k^2 * n_l
    // e.g. n_(0,0,0) * n_(1,0,1) * n_(0,1,1) * n_(1,1,0)
    double m6; // sum_NNN\, n_i^2 * n_i+2 (i=x,y,z)

    double x_1; // sum\, n_i^2*v_i
    double x0; // sum_NN\, n_i^2*v_i * n_i+1 (i=x,y,z)
    double x1; // sum_{diagonal in plane}\, n_i^2*v_i * n_j
    // e.g. n_(0,0,0) * n_(1,1,0)
    double x2; // sum_{cubic diagonal}\, n_i^2*v_i * n_j
    // e.g. n_(0,0,0) * n_(1,1,1)
    double x3; // sum_{principal planes}\, n_i^2*v_i * n_j * n_k^2*v_k * n_l
    // e.g. n_(0,0,0) * n_(0,1,0) * n_(1,0,0) * n_(1,1,0)
    double x4; // sum_{diagonal planes}\, n_i^2*v_i * n_j * n_k^2*v_k * n_l
    // e.g. n_(0,0,0) * n_(1,0,0) * n_(0,1,1) * n_(1,1,1)
    double x5; // sum_{tetrahedral vertices}\, n_i^2*v_i * n_j * n_k^2*v_k * n_l
    // e.g. n_(0,0,0) * n_(1,0,1) * n_(0,1,1) * n_(1,1,0)
    double x6; // sum_NNN\, n_i^2*v_i * n_i+2 (i=x,y,z)

    double c_1; // sum\, n_i*v_i*c_i
    double c0; // sum_NN\, n_i*v_i*c_i * n_i+1 (i=x,y,z)
    double c1; // sum_{diagonal in plane}\, n_i*v_i*c_i * n_j
    // e.g. n_(0,0,0) * n_(1,1,0)
    double c2; // sum_{cubic diagonal}\, n_i*v_i*c_i * n_j
    // e.g. n_(0,0,0) * n_(1,1,1)
    double c3; // sum_{principal planes}\, n_i*v_i*c_i * n_j * n_k*v_k*c_k * n_l
    // e.g. n_(0,0,0) * n_(0,1,0) * n_(1,0,0) * n_(1,1,0)
    double c4; // sum_{diagonal planes}\, n_i*v_i*c_i * n_j * n_k*v_k*c_k * n_l
    // e.g. n_(0,0,0) * n_(1,0,0) * n_(0,1,1) * n_(1,1,1)
    double c5; // sum_{tetrahedral vertices}\, n_i*v_i*c_i * n_j * n_k*v_k*c_k * n_l
    // e.g. n_(0,0,0) * n_(1,0,1) * n_(0,1,1) * n_(1,1,0)
    double c6; // sum_NNN\, n_i*v_i*c_i * n_i+2 (i=x,y,z)

    double n_1vec[runpts];
    double n0vec[runpts];
    double n1vec[runpts];
    double n2vec[runpts];
    double n3vec[runpts];
    double n4vec[runpts];
    double n5vec[runpts];
    double n6vec[runpts];

    double u_1vec[runpts];
    double u0vec[runpts];
    double u1vec[runpts];
    double u2vec[runpts];
    double u3vec[runpts];
    double u4vec[runpts];
    double u5vec[runpts];
    double u6vec[runpts];

    double m_1vec[runpts];
    double m0vec[runpts];
    double m1vec[runpts];
    double m2vec[runpts];
    double m3vec[runpts];
    double m4vec[runpts];
    double m5vec[runpts];
    double m6vec[runpts];

    double x_1vec[runpts];
    double x0vec[runpts];
    double x1vec[runpts];
    double x2vec[runpts];
    double x3vec[runpts];
    double x4vec[runpts];
    double x5vec[runpts];
    double x6vec[runpts];

    double c_1vec[runpts];
    double c0vec[runpts];
    double c1vec[runpts];
    double c2vec[runpts];
    double c3vec[runpts];
    double c4vec[runpts];
    double c5vec[runpts];
    double c6vec[runpts];

    int celli,cellj,cellk,celll;

    double nn;
    double vv;
    double cc;

    double nn_i,nn_j,nn_k,nn_l;
    double vv_i,vv_j,vv_k,vv_l;
    double cc_i,cc_j,cc_k,cc_l;

    for (i=0; i<(int)(pow(L/b,dim)); i++)
    {
        ni[i] = 0.0;
        vi[i] = 0.0;
        ci[i] = 0.0;
    }

    int runpt;

    for (runpt=0; runpt<runpts; runpt++)
    {
        if (b>1)
        {
            bshift = rand() % b;
        }
        else
        {
            bshift = 0;
        }

        // read in data

        filenameStream << "./netlib/n" << runpt << ".dat";
        filename = filenameStream.str();

        fin.open(filename.c_str());

        j = 0;
        while( getline(fin,line) )
        {
            n[j] = atoi (line.c_str());
            j++;
        }

        fin.close();
        filenameStream.str("");

        filenameStream << "./netlib/via" << runpt << ".dat";
        filename = filenameStream.str();

        fin.open(filename.c_str());

        j = 0;
        while( getline(fin,line) )
        {
            via[j] = atof (line.c_str());
            j++;
        }

        fin.close();
        filenameStream.str("");

        filenameStream << "./netlib/c" << runpt << ".dat";
        filename = filenameStream.str();

        fin.open(filename.c_str());

        j = 0;
        while( getline(fin,line) )
        {
            c[j] = atof (line.c_str());
            j++;
        }

        fin.close();
        filenameStream.str("");

        for (i=0; i<nL; i++)
        {
            v[i] = 0.0;
            for (j=0; j<nLfine; j++)
            {
                v[i] = v[i]+via[i*nLfine+j];
            }
        }

        i = 0; j = 0; k = 0;

        ib = 0; jb = 0; kb = 0;

        for (i=0; i<(int)(pow(L/b,dim)); i++)
        {
            ni[i] = 0.0;
            vi[i] = 0.0;
            ci[i] = 0.0;
        }

        n_1 = 0.0; n0 = 0.0; n1 = 0.0; n2 = 0.0; n3 = 0.0; n4 = 0.0; n5 = 0.0; n6 = 0.0;
        u_1 = 0.0; u0 = 0.0; u1 = 0.0; u2 = 0.0; u3 = 0.0; u4 = 0.0; u5 = 0.0; u6 = 0.0;

        m_1 = 0.0; m0 = 0.0; m1 = 0.0; m2 = 0.0; m3 = 0.0; m4 = 0.0; m5 = 0.0; m6 = 0.0;
        x_1 = 0.0; x0 = 0.0; x1 = 0.0; x2 = 0.0; x3 = 0.0; x4 = 0.0; x5 = 0.0; x6 = 0.0;

        c_1 = 0.0; c0 = 0.0; c1 = 0.0; c2 = 0.0; c3 = 0.0; c4 = 0.0; c5 = 0.0; c6 = 0.0;

        nn = 0.0;
        vv = 0.0;
        cc = 0.0;

        nn_i = 0.0; nn_j = 0.0; nn_k = 0.0; nn_l = 0.0;
        vv_i = 0.0; vv_j = 0.0; vv_k = 0.0; vv_l = 0.0;
        cc_i = 0.0; cc_j = 0.0; cc_k = 0.0; cc_l = 0.0;

        celli = 0; cellj = 0; cellk = 0; celll = 0;

    #pragma omp parallel for shared (n,v,c,ni,vi,ci) private (i,j,k,celli,nn,vv,cc,ib,jb,kb)
        for (i=0; i<(int)(L/b); i++)
        {
            for (j=0; j<(int)(L/b); j++)
            {
                for (k=0; k<(int)(L/b); k++)
                {
            // calculate average n for this cell block
            // calculate average excluded volume for this cell block

                    celli = (int)(i*pow(L/b,2)+j*(L/b)+k);

        #pragma omp atomic write
                    ni[celli] = 0.0;
        #pragma omp atomic write
                    vi[celli] = 0.0;
        #pragma omp atomic write
                    ci[celli] = 0.0;

                    nn = 0.0;
                    vv = 0.0;
                    cc = 0.0;

                    for (ib=bshift; ib<(bshift+b); ib++)
                    {
                        for (jb=bshift; jb<(bshift+b); jb++)
                        {
                            for (kb=bshift; kb<(bshift+b); kb++)
                            {

                                if ((i*b+ib)>(L-1))
                                {
                                    if ((j*b+jb)>(L-1))
                                    {
                                        if ((k*b+kb)>(L-1))
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb-L)*L+k*b+kb-L)];
        #pragma omp atomic read
                                            vv = v[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb-L)*L+k*b+kb-L)];
        #pragma omp atomic read
                                            cc = c[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb-L)*L+k*b+kb-L)];
                                        }
                                        else
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb-L)*L+k*b+kb)];
        #pragma omp atomic read
                                            vv = v[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb-L)*L+k*b+kb)];
        #pragma omp atomic read
                                            cc = c[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb-L)*L+k*b+kb)];
                                        }
                                    }
                                    else
                                    {
                                        if ((k*b+kb)>(L-1))
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb)*L+k*b+kb-L)];
        #pragma omp atomic read
                                            vv = v[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb)*L+k*b+kb-L)];
        #pragma omp atomic read
                                            cc = c[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb)*L+k*b+kb-L)];
                                        }
                                        else
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb)*L+k*b+kb)];
        #pragma omp atomic read
                                            vv = v[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb)*L+k*b+kb)];
        #pragma omp atomic read
                                            cc = c[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb)*L+k*b+kb)];
                                        }

                                    }
                                }
                                else
                                {
                                    if ((j*b+jb)>(L-1))
                                    {
                                        if ((k*b+kb)>(L-1))
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib)*pow(L,2)+(j*b+jb-L)*L+k*b+kb-L)];
        #pragma omp atomic read
                                            vv = v[(int)((i*b+ib)*pow(L,2)+(j*b+jb-L)*L+k*b+kb-L)];
        #pragma omp atomic read
                                            cc = c[(int)((i*b+ib)*pow(L,2)+(j*b+jb-L)*L+k*b+kb-L)];
                                        }
                                        else
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib)*pow(L,2)+(j*b+jb-L)*L+k*b+kb)];
        #pragma omp atomic read
                                            vv = v[(int)((i*b+ib)*pow(L,2)+(j*b+jb-L)*L+k*b+kb)];
        #pragma omp atomic read
                                            cc = c[(int)((i*b+ib)*pow(L,2)+(j*b+jb-L)*L+k*b+kb)];

                                        }
                                    }
                                    else
                                    {
                                        if ((k*b+kb)>(L-1))
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib)*pow(L,2)+(j*b+jb)*L+k*b+kb-L)];
        #pragma omp atomic read
                                            vv = v[(int)((i*b+ib)*pow(L,2)+(j*b+jb)*L+k*b+kb-L)];
        #pragma omp atomic read
                                            cc = c[(int)((i*b+ib)*pow(L,2)+(j*b+jb)*L+k*b+kb-L)];
                                        }
                                        else
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib)*pow(L,2)+(j*b+jb)*L+k*b+kb)];
        #pragma omp atomic read
                                            vv = v[(int)((i*b+ib)*pow(L,2)+(j*b+jb)*L+k*b+kb)];
        #pragma omp atomic read
                                            cc = c[(int)((i*b+ib)*pow(L,2)+(j*b+jb)*L+k*b+kb)];
                                        }

                                    }

                                }

        #pragma omp atomic write
                                ni[celli] = ni[celli] + nn;
        #pragma omp atomic write
                                vi[celli] = vi[celli] + vv;
        #pragma omp atomic write
                                ci[celli] = ci[celli] + cc;

                            }
                        }
                    }
                }
            }
        }

    #pragma omp parallel for shared (n,v,c,ni,vi,ci) private (i,j,k,celli,nn,vv,cc) reduction(+:n_1,u_1,m_1,x_1,c_1)
        for (i=0; i<(int)(L/b); i++)
        {
            for (j=0; j<(int)(L/b); j++)
            {
                for (k=0; k<(int)(L/b); k++)
                {
                    celli = (int)(i*pow(L/b,2)+j*(L/b)+k);

        #pragma omp atomic read
                    nn = ni[celli];
        #pragma omp atomic read
                    vv = vi[celli];
        #pragma omp atomic read
                    cc = ci[celli];

        #pragma omp atomic write
                    ni[celli] = nn/(double)(pow(b,dim));
        #pragma omp atomic write
                    vi[celli] = vv/(double)(pow(b,dim));
        #pragma omp atomic write
                    ci[celli] = cc/(double)(pow(b,dim));

                    n_1 += nn/(double)(pow(b,dim));
                    u_1 += nn*vv/(double)(pow(b,2*dim));
                    m_1 += nn*nn/(double)(pow(b,2*dim));
                    x_1 += nn*nn*vv/(double)(pow(b,3*dim));
                    c_1 += nn*vv*cc/(double)(pow(b,3*dim));

                }
            }
        }

        ndump.open("ni.dat", std::ios_base::app);

        for (i=0; i<(int)(pow(L/b,dim)); i++)
        {
            ndump << ni[i] << "\n";
        }
        ndump.close();

        ndump.open("vi.dat", std::ios_base::app);

        for (i=0; i<(int)(pow(L/b,dim)); i++)
        {
            ndump << vi[i] << "\n";
        }
        ndump.close();

        ndump.open("ci.dat", std::ios_base::app);

        for (i=0; i<(int)(pow(L/b,dim)); i++)
        {
            ndump << ci[i] << "\n";
        }
        ndump.close();

    #pragma omp parallel for shared (ni,vi,ci) private (i,j,k,celli,cellj,nn_i,nn_j,vv_i,vv_j,cc_i,cc_j) reduction(+:n0,u0,m0,x0,c0)
        for (i=0; i<(int)(L/b); i++)
        {
            for (j=0; j<(int)(L/b); j++)
            {
                for (k=0; k<(int)(L/b); k++)
                {
                    // sum_NN\, n_i * n_i+1 (i=x,y,z)

                    celli = (int)(i*pow(L/b,2)+j*(L/b)+k);
        #pragma omp atomic read
                    nn_i = ni[celli];
        #pragma omp atomic read
                    vv_i = vi[celli];
        #pragma omp atomic read
                    cc_i = ci[celli];

                    cellj = get_cell(i,j,k,1,0,0,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n0 += nn_i*nn_j;
                    u0 += nn_i*vv_i*nn_j;
                    m0 += nn_i*nn_i*nn_j;
                    x0 += nn_i*nn_i*vv_i*nn_j;
                    c0 += nn_i*vv_i*cc_i*nn_j;

                    cellj = get_cell(i,j,k,0,1,0,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n0 += nn_i*nn_j;
                    u0 += nn_i*vv_i*nn_j;
                    m0 += nn_i*nn_i*nn_j;
                    x0 += nn_i*nn_i*vv_i*nn_j;
                    c0 += nn_i*vv_i*cc_i*nn_j;

                    cellj = get_cell(i,j,k,0,0,1,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n0 += nn_i*nn_j;
                    u0 += nn_i*vv_i*nn_j;
                    m0 += nn_i*nn_i*nn_j;
                    x0 += nn_i*nn_i*vv_i*nn_j;
                    c0 += nn_i*vv_i*cc_i*nn_j;
                }
            }
        }

    #pragma omp parallel for shared (ni,vi,ci) private (i,j,k,celli,cellj,nn_i,nn_j,vv_i,vv_j,cc_i,cc_j) reduction(+:n1,u1,m1,x1,c1)
        for (i=0; i<(int)(L/b); i++)
        {
            for (j=0; j<(int)(L/b); j++)
            {
                for (k=0; k<(int)(L/b); k++)
                {
                    // sum_{diagonal in plane}\, n_i * n_j
                    // e.g. n_(0,0,0) * n_(1,1,0)

                    celli = (int)(i*pow(L/b,2)+j*(L/b)+k);
        #pragma omp atomic read
                    nn_i = ni[celli];
        #pragma omp atomic read
                    vv_i = vi[celli];
        #pragma omp atomic read
                    cc_i = ci[celli];

                    cellj = get_cell(i,j,k,1,1,0,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n1 += nn_i*nn_j; 
                    u1 += nn_i*vv_i*nn_j;
                    m1 += nn_i*nn_i*nn_j;
                    x1 += nn_i*nn_i*vv_i*nn_j;
                    c1 += nn_i*vv_i*cc_i*nn_j;

                    cellj = get_cell(i,j,k,1,-1,0,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n1 += nn_i*nn_j; 
                    u1 += nn_i*vv_i*nn_j;
                    m1 += nn_i*nn_i*nn_j;
                    x1 += nn_i*nn_i*vv_i*nn_j;
                    c1 += nn_i*vv_i*cc_i*nn_j;

                    cellj = get_cell(i,j,k,1,0,1,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n1 += nn_i*nn_j; 
                    u1 += nn_i*vv_i*nn_j;
                    m1 += nn_i*nn_i*nn_j;
                    x1 += nn_i*nn_i*vv_i*nn_j;
                    c1 += nn_i*vv_i*cc_i*nn_j;

                    cellj = get_cell(i,j,k,1,0,-1,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n1 += nn_i*nn_j; 
                    u1 += nn_i*vv_i*nn_j;
                    m1 += nn_i*nn_i*nn_j;
                    x1 += nn_i*nn_i*vv_i*nn_j;
                    c1 += nn_i*vv_i*cc_i*nn_j;

                    cellj = get_cell(i,j,k,0,1,1,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n1 += nn_i*nn_j; 
                    u1 += nn_i*vv_i*nn_j;
                    m1 += nn_i*nn_i*nn_j;
                    x1 += nn_i*nn_i*vv_i*nn_j;
                    c1 += nn_i*vv_i*cc_i*nn_j;

                    cellj = get_cell(i,j,k,0,1,-1,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n1 += nn_i*nn_j; 
                    u1 += nn_i*vv_i*nn_j;
                    m1 += nn_i*nn_i*nn_j;
                    x1 += nn_i*nn_i*vv_i*nn_j;
                    c1 += nn_i*vv_i*cc_i*nn_j;
                }
            }
        }

    #pragma omp parallel for shared (ni,vi,ci) private (i,j,k,celli,cellj,nn_i,nn_j,vv_i,vv_j,cc_i,cc_j) reduction(+:n2,u2,m2,x2,c2)
        for (i=0; i<(int)(L/b); i++)
        {
            for (j=0; j<(int)(L/b); j++)
            {
                for (k=0; k<(int)(L/b); k++)
                {
                    // sum_{cubic diagonal}\, n_i * n_j
                    // e.g. n_(0,0,0) * n_(1,1,1)

                    celli = (int)(i*pow(L/b,2)+j*(L/b)+k);
        #pragma omp atomic read
                    nn_i = ni[celli];
        #pragma omp atomic read
                    vv_i = vi[celli];
        #pragma omp atomic read
                    cc_i = ci[celli];

                    cellj = get_cell(i,j,k,1,1,1,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n2 += nn_i*nn_j; 
                    u2 += nn_i*vv_i*nn_j;
                    m2 += nn_i*nn_i*nn_j;
                    x2 += nn_i*nn_i*vv_i*nn_j;
                    c2 += nn_i*vv_i*cc_i*nn_j;

                    cellj = get_cell(i,j,k,1,1,-1,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n2 += nn_i*nn_j; 
                    u2 += nn_i*vv_i*nn_j;
                    m2 += nn_i*nn_i*nn_j;
                    x2 += nn_i*nn_i*vv_i*nn_j;
                    c2 += nn_i*vv_i*cc_i*nn_j;

                    cellj = get_cell(i,j,k,1,-1,-1,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n2 += nn_i*nn_j; 
                    u2 += nn_i*vv_i*nn_j;
                    m2 += nn_i*nn_i*nn_j;
                    x2 += nn_i*nn_i*vv_i*nn_j;
                    c2 += nn_i*vv_i*cc_i*nn_j;

                    cellj = get_cell(i,j,k,1,-1,1,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n2 += nn_i*nn_j; 
                    u2 += nn_i*vv_i*nn_j;
                    m2 += nn_i*nn_i*nn_j;
                    x2 += nn_i*nn_i*vv_i*nn_j;
                    c2 += nn_i*vv_i*cc_i*nn_j;
                }
            }
        }

    #pragma omp parallel for shared (ni,vi,ci) private (i,j,k,celli,cellj,cellk,celll,nn,nn_i,nn_j,nn_k,nn_l,vv_i,vv_j,vv_k,vv_l,cc_i,cc_j,cc_k,cc_l) reduction(+:n3,u3,m3,x3,c3)
        for (i=0; i<(int)(L/b); i++)
        {
            for (j=0; j<(int)(L/b); j++)
            {
                for (k=0; k<(int)(L/b); k++)
                {
                    // sum_{principal planes}\, n_i * n_j * n_k * n_l
                    // e.g. n_(0,0,0) * n_(0,1,0) * n_(1,0,0) * n_(1,1,0)

                    celli = (int)(i*pow(L/b,2)+j*(L/b)+k);
        #pragma omp atomic read
                    nn_i = ni[celli];
        #pragma omp atomic read
                    vv_i = vi[celli];
        #pragma omp atomic read
                    cc_i = ci[celli];

                    cellj = get_cell(i,j,k,1,0,0,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    cellk = get_cell(i,j,k,0,1,0,b);

        #pragma omp atomic read
                    nn_k = ni[cellk];
        #pragma omp atomic read
                    vv_k = vi[cellk];
        #pragma omp atomic read
                    cc_k = ci[cellk];

                    celll = get_cell(i,j,k,1,1,0,b);

        #pragma omp atomic read
                    nn_l = ni[celll];
        #pragma omp atomic read
                    vv_l = vi[celll];
        #pragma omp atomic read
                    cc_l = ci[celll];

                    n3 += nn_i*nn_j*nn_k*nn_l; 
                    u3 += nn_i*vv_i*nn_j*nn_k*vv_k*nn_l;
                    m3 += nn_i*nn_i*nn_j*nn_k*nn_k*nn_l;
                    x3 += nn_i*nn_i*vv_i*nn_j*nn_k*nn_k*vv_k*nn_l;
                    c3 += nn_i*vv_i*cc_i*nn_j*nn_k*vv_k*cc_k*nn_l;

                    cellj = get_cell(i,j,k,1,0,0,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    cellk = get_cell(i,j,k,0,0,1,b);

        #pragma omp atomic read
                    nn_k = ni[cellk];
        #pragma omp atomic read
                    vv_k = vi[cellk];
        #pragma omp atomic read
                    cc_k = ci[cellk];

                    celll = get_cell(i,j,k,1,0,1,b);

        #pragma omp atomic read
                    nn_l = ni[celll];
        #pragma omp atomic read
                    vv_l = vi[celll];
        #pragma omp atomic read
                    cc_l = ci[celll];

                    n3 += nn_i*nn_j*nn_k*nn_l; 
                    u3 += nn_i*vv_i*nn_j*nn_k*vv_k*nn_l;
                    m3 += nn_i*nn_i*nn_j*nn_k*nn_k*nn_l;
                    x3 += nn_i*nn_i*vv_i*nn_j*nn_k*nn_k*vv_k*nn_l;
                    c3 += nn_i*vv_i*cc_i*nn_j*nn_k*vv_k*cc_k*nn_l;

                    cellj = get_cell(i,j,k,0,1,0,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    cellk = get_cell(i,j,k,0,0,1,b);

        #pragma omp atomic read
                    nn_k = ni[cellk];
        #pragma omp atomic read
                    vv_k = vi[cellk];
        #pragma omp atomic read
                    cc_k = ci[cellk];

                    celll = get_cell(i,j,k,0,1,1,b);

        #pragma omp atomic read
                    nn_l = ni[celll];
        #pragma omp atomic read
                    vv_l = vi[celll];
        #pragma omp atomic read
                    cc_l = ci[celll];

                    n3 += nn_i*nn_j*nn_k*nn_l; 
                    u3 += nn_i*vv_i*nn_j*nn_k*vv_k*nn_l;
                    m3 += nn_i*nn_i*nn_j*nn_k*nn_k*nn_l;
                    x3 += nn_i*nn_i*vv_i*nn_j*nn_k*nn_k*vv_k*nn_l;
                    c3 += nn_i*vv_i*cc_i*nn_j*nn_k*vv_k*cc_k*nn_l;
                }
            }
        }

    #pragma omp parallel for shared (ni,vi,ci) private (i,j,k,celli,cellj,cellk,celll,nn,nn_i,nn_j,nn_k,nn_l,vv_i,vv_j,vv_k,vv_l,cc_i,cc_j,cc_k,cc_l) reduction(+:n4,u4,m4,x4,c4)
        for (i=0; i<(int)(L/b); i++)
        {
            for (j=0; j<(int)(L/b); j++)
            {
                for (k=0; k<(int)(L/b); k++)
                {
                    // sum_{diagonal planes}\, n_i * n_j * n_k * n_l
                    // e.g. n_(0,0,0) * n_(1,0,0) * n_(0,1,1) * n_(1,1,1)

                    celli = (int)(i*pow(L/b,2)+j*(L/b)+k);
        #pragma omp atomic read
                    nn_i = ni[celli];
        #pragma omp atomic read
                    vv_i = vi[celli];
        #pragma omp atomic read
                    cc_i = ci[celli];

                    cellj = get_cell(i,j,k,1,0,0,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    cellk = get_cell(i,j,k,0,1,1,b);

        #pragma omp atomic read
                    nn_k = ni[cellk];
        #pragma omp atomic read
                    vv_k = vi[cellk];
        #pragma omp atomic read
                    cc_k = ci[cellk];

                    celll = get_cell(i,j,k,1,1,1,b);

        #pragma omp atomic read
                    nn_l = ni[celll];
        #pragma omp atomic read
                    vv_l = vi[celll];
        #pragma omp atomic read
                    cc_l = ci[celll];

                    n4 += nn_i*nn_j*nn_k*nn_l; 
                    u4 += nn_i*vv_i*nn_j*nn_k*vv_k*nn_l;
                    m4 += nn_i*nn_i*nn_j*nn_k*nn_k*nn_l;
                    x4 += nn_i*nn_i*vv_i*nn_j*nn_k*nn_k*vv_k*nn_l;
                    c4 += nn_i*vv_i*cc_i*nn_j*nn_k*vv_k*cc_k*nn_l;

                    cellj = get_cell(i,j,k,1,0,0,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    cellk = get_cell(i,j,k,0,1,-1,b);

        #pragma omp atomic read
                    nn_k = ni[cellk];
        #pragma omp atomic read
                    vv_k = vi[cellk];
        #pragma omp atomic read
                    cc_k = ci[cellk];

                    celll = get_cell(i,j,k,1,1,-1,b);

        #pragma omp atomic read
                    nn_l = ni[celll];
        #pragma omp atomic read
                    vv_l = vi[celll];
        #pragma omp atomic read
                    cc_l = ci[celll];

                    n4 += nn_i*nn_j*nn_k*nn_l; 
                    u4 += nn_i*vv_i*nn_j*nn_k*vv_k*nn_l;
                    m4 += nn_i*nn_i*nn_j*nn_k*nn_k*nn_l;
                    x4 += nn_i*nn_i*vv_i*nn_j*nn_k*nn_k*vv_k*nn_l;
                    c4 += nn_i*vv_i*cc_i*nn_j*nn_k*vv_k*cc_k*nn_l;

                    cellj = get_cell(i,j,k,0,1,0,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    cellk = get_cell(i,j,k,1,0,1,b);

        #pragma omp atomic read
                    nn_k = ni[cellk];
        #pragma omp atomic read
                    vv_k = vi[cellk];
        #pragma omp atomic read
                    cc_k = ci[cellk];

                    celll = get_cell(i,j,k,1,1,1,b);

        #pragma omp atomic read
                    nn_l = ni[celll];
        #pragma omp atomic read
                    vv_l = vi[celll];
        #pragma omp atomic read
                    cc_l = ci[celll];

                    n4 += nn_i*nn_j*nn_k*nn_l; 
                    u4 += nn_i*vv_i*nn_j*nn_k*vv_k*nn_l;
                    m4 += nn_i*nn_i*nn_j*nn_k*nn_k*nn_l;
                    x4 += nn_i*nn_i*vv_i*nn_j*nn_k*nn_k*vv_k*nn_l;
                    c4 += nn_i*vv_i*cc_i*nn_j*nn_k*vv_k*cc_k*nn_l;

                    cellj = get_cell(i,j,k,0,1,0,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    cellk = get_cell(i,j,k,1,0,-1,b);

        #pragma omp atomic read
                    nn_k = ni[cellk];
        #pragma omp atomic read
                    vv_k = vi[cellk];
        #pragma omp atomic read
                    cc_k = ci[cellk];

                    celll = get_cell(i,j,k,1,1,-1,b);

        #pragma omp atomic read
                    nn_l = ni[celll];
        #pragma omp atomic read
                    vv_l = vi[celll];
        #pragma omp atomic read
                    cc_l = ci[celll];

                    n4 += nn_i*nn_j*nn_k*nn_l; 
                    u4 += nn_i*vv_i*nn_j*nn_k*vv_k*nn_l;
                    m4 += nn_i*nn_i*nn_j*nn_k*nn_k*nn_l;
                    x4 += nn_i*nn_i*vv_i*nn_j*nn_k*nn_k*vv_k*nn_l;
                    c4 += nn_i*vv_i*cc_i*nn_j*nn_k*vv_k*cc_k*nn_l;

                    cellj = get_cell(i,j,k,0,0,1,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    cellk = get_cell(i,j,k,1,1,0,b);

        #pragma omp atomic read
                    nn_k = ni[cellk];
        #pragma omp atomic read
                    vv_k = vi[cellk];
        #pragma omp atomic read
                    cc_k = ci[cellk];

                    celll = get_cell(i,j,k,1,1,1,b);

        #pragma omp atomic read
                    nn_l = ni[celll];
        #pragma omp atomic read
                    vv_l = vi[celll];
        #pragma omp atomic read
                    cc_l = ci[celll];

                    n4 += nn_i*nn_j*nn_k*nn_l; 
                    u4 += nn_i*vv_i*nn_j*nn_k*vv_k*nn_l;
                    m4 += nn_i*nn_i*nn_j*nn_k*nn_k*nn_l;
                    x4 += nn_i*nn_i*vv_i*nn_j*nn_k*nn_k*vv_k*nn_l;
                    c4 += nn_i*vv_i*cc_i*nn_j*nn_k*vv_k*cc_k*nn_l;

                    cellj = get_cell(i,j,k,0,0,1,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    cellk = get_cell(i,j,k,1,-1,0,b);

        #pragma omp atomic read
                    nn_k = ni[cellk];
        #pragma omp atomic read
                    vv_k = vi[cellk];
        #pragma omp atomic read
                    cc_k = ci[cellk];

                    celll = get_cell(i,j,k,1,-1,1,b);

        #pragma omp atomic read
                    nn_l = ni[celll];
        #pragma omp atomic read
                    vv_l = vi[celll];
        #pragma omp atomic read
                    cc_l = ci[celll];

                    n4 += nn_i*nn_j*nn_k*nn_l; 
                    u4 += nn_i*vv_i*nn_j*nn_k*vv_k*nn_l;
                    m4 += nn_i*nn_i*nn_j*nn_k*nn_k*nn_l;
                    x4 += nn_i*nn_i*vv_i*nn_j*nn_k*nn_k*vv_k*nn_l;
                    c4 += nn_i*vv_i*cc_i*nn_j*nn_k*vv_k*cc_k*nn_l;
                }
            }
        }

    #pragma omp parallel for shared (ni,vi,ci) private (i,j,k,celli,cellj,cellk,celll,nn,nn_i,nn_j,nn_k,nn_l,vv_i,vv_j,vv_k,vv_l,cc_i,cc_j,cc_k,cc_l) reduction(+:n5,u5,m5,x5,c5)
        for (i=0; i<(int)(L/b); i++)
        {
            for (j=0; j<(int)(L/b); j++)
            {
                for (k=0; k<(int)(L/b); k++)
                {
                    // sum_{tetrahedral vertices}\, n_i * n_j * n_k * n_l
                    // e.g. n_(0,0,0) * n_(1,0,1) * n_(0,1,1) * n_(1,1,0)

                    celli = (int)(i*pow(L/b,2)+j*(L/b)+k);
        #pragma omp atomic read
                    nn_i = ni[celli];
        #pragma omp atomic read
                    vv_i = vi[celli];
        #pragma omp atomic read
                    cc_i = ci[celli];

                    cellj = get_cell(i,j,k,1,0,1,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    cellk = get_cell(i,j,k,0,1,1,b);

        #pragma omp atomic read
                    nn_k = ni[cellk];
        #pragma omp atomic read
                    vv_k = vi[cellk];
        #pragma omp atomic read
                    cc_k = ci[cellk];

                    celll = get_cell(i,j,k,1,1,0,b);

        #pragma omp atomic read
                    nn_l = ni[celll];
        #pragma omp atomic read
                    vv_l = vi[celll];
        #pragma omp atomic read
                    cc_l = ci[celll];

                    n5 += nn_i*nn_j*nn_k*nn_l; 
                    u5 += nn_i*vv_i*nn_j*nn_k*vv_k*nn_l;
                    m5 += nn_i*nn_i*nn_j*nn_k*nn_k*nn_l;
                    x5 += nn_i*nn_i*vv_i*nn_j*nn_k*nn_k*vv_k*nn_l;
                    c5 += nn_i*vv_i*cc_i*nn_j*nn_k*vv_k*cc_k*nn_l;

                    cellj = get_cell(i,j,k,1,0,1,b);

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    cellk = get_cell(i,j,k,0,-1,1,b);

        #pragma omp atomic read
                    nn_k = ni[cellk];
        #pragma omp atomic read
                    vv_k = vi[cellk];
        #pragma omp atomic read
                    cc_k = ci[cellk];

                    celll = get_cell(i,j,k,1,-1,0,b);

        #pragma omp atomic read
                    nn_l = ni[celll];
        #pragma omp atomic read
                    vv_l = vi[celll];
        #pragma omp atomic read
                    cc_l = ci[celll];

                    n5 += nn_i*nn_j*nn_k*nn_l; 
                    u5 += nn_i*vv_i*nn_j*nn_k*vv_k*nn_l;
                    m5 += nn_i*nn_i*nn_j*nn_k*nn_k*nn_l;
                    x5 += nn_i*nn_i*vv_i*nn_j*nn_k*nn_k*vv_k*nn_l;
                    c5 += nn_i*vv_i*cc_i*nn_j*nn_k*vv_k*cc_k*nn_l;
                }
            }
        }

    #pragma omp parallel for shared (ni,vi,ci) private (i,j,k,celli,cellj,nn_i,nn_j,vv_i,vv_j,cc_i,cc_j) reduction(+:n6,u6,m6,x6,c6)
        for (i=0; i<(int)(L/b); i++)
        {
            for (j=0; j<(int)(L/b); j++)
            {
                for (k=0; k<(int)(L/b); k++)
                {
                    // sum_NNN\, n_i * n_i+2 (i=x,y,z)

                    celli = (int)(i*pow(L/b,2)+j*(L/b)+k);
        #pragma omp atomic read
                    nn_i = ni[celli];
        #pragma omp atomic read
                    vv_i = vi[celli];
        #pragma omp atomic read
                    cc_i = ci[celli];

                    //cellj = get_cell(i,j,k,2,0,0,b); 
                    if (i==((L/b)-1))
                    {
                        cellj = get_cell(0,j,k,1,0,0,b);
                    }
                    else
                    {
                        cellj = get_cell(i+1,j,k,1,0,0,b);
                    }

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n6 += nn_i*nn_j; 
                    u6 += nn_i*vv_i*nn_j;
                    m6 += nn_i*nn_i*nn_j;
                    x6 += nn_i*nn_i*vv_i*nn_j;
                    c6 += nn_i*vv_i*cc_i*nn_j;

                    //cellj = get_cell(i,j,k,0,2,0,b); 
                    if (j==((L/b)-1))
                    {
                        cellj = get_cell(i,0,k,0,1,0,b);
                    }
                    else
                    {
                        cellj = get_cell(i,j+1,k,0,1,0,b);
                    }

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n6 += nn_i*nn_j; 
                    u6 += nn_i*vv_i*nn_j;
                    m6 += nn_i*nn_i*nn_j;
                    x6 += nn_i*nn_i*vv_i*nn_j;
                    c6 += nn_i*vv_i*cc_i*nn_j;

                    //cellj = get_cell(i,j,k,0,0,2,b); 
                    if (k==((L/b)-1))
                    {
                        cellj = get_cell(i,j,0,0,0,1,b);
                    }
                    else
                    {
                        cellj = get_cell(i,j,k+1,0,0,1,b);
                    }

        #pragma omp atomic read
                    nn_j = ni[cellj];
        #pragma omp atomic read
                    vv_j = vi[cellj];
        #pragma omp atomic read
                    cc_j = ci[cellj];

                    n6 += nn_i*nn_j; 
                    u6 += nn_i*vv_i*nn_j;
                    m6 += nn_i*nn_i*nn_j;
                    x6 += nn_i*nn_i*vv_i*nn_j;
                    c6 += nn_i*vv_i*cc_i*nn_j;
                }
            }
        }

        n_1vec[runpt] = n_1;
        n0vec[runpt] = n0;
        n1vec[runpt] = n1;
        n2vec[runpt] = n2;
        n3vec[runpt] = n3;
        n4vec[runpt] = n4;
        n5vec[runpt] = n5;
        n6vec[runpt] = n6;

        u_1vec[runpt] = u_1;
        u0vec[runpt] = u0;
        u1vec[runpt] = u1;
        u2vec[runpt] = u2;
        u3vec[runpt] = u3;
        u4vec[runpt] = u4;
        u5vec[runpt] = u5;
        u6vec[runpt] = u6;

        m_1vec[runpt] = m_1;
        m0vec[runpt] = m0;
        m1vec[runpt] = m1;
        m2vec[runpt] = m2;
        m3vec[runpt] = m3;
        m4vec[runpt] = m4;
        m5vec[runpt] = m5;
        m6vec[runpt] = m6;

        x_1vec[runpt] = x_1;
        x0vec[runpt] = x0;
        x1vec[runpt] = x1;
        x2vec[runpt] = x2;
        x3vec[runpt] = x3;
        x4vec[runpt] = x4;
        x5vec[runpt] = x5;
        x6vec[runpt] = x6;

        c_1vec[runpt] = c_1;
        c0vec[runpt] = c0;
        c1vec[runpt] = c1;
        c2vec[runpt] = c2;
        c3vec[runpt] = c3;
        c4vec[runpt] = c4;
        c5vec[runpt] = c5;
        c6vec[runpt] = c6;
    }

    // save results to disk

    ndump.open("nvec.dat", std::ios_base::app);

    for (runpt=0; runpt<runpts; runpt++)
    {
        ndump << n_1vec[runpt] << "\n" << n0vec[runpt] << "\n" << n1vec[runpt] << "\n" << 
                 n2vec[runpt] << "\n" << n3vec[runpt] << "\n" <<
                 n4vec[runpt] << "\n" << n5vec[runpt] << "\n" <<
                 n6vec[runpt] << "\n";
    }

    ndump.close();

    ndump.open("uvec.dat", std::ios_base::app);

    for (runpt=0; runpt<runpts; runpt++)
    {
        ndump << u_1vec[runpt] << "\n" << u0vec[runpt] << "\n" << u1vec[runpt] << "\n" << 
                 u2vec[runpt] << "\n" << u3vec[runpt] << "\n" <<
                 u4vec[runpt] << "\n" << u5vec[runpt] << "\n" <<
                 u6vec[runpt] << "\n";
    }

    ndump.close();

    ndump.open("mvec.dat", std::ios_base::app);

    for (runpt=0; runpt<runpts; runpt++)
    {
        ndump << m_1vec[runpt] << "\n" << m0vec[runpt] << "\n" << m1vec[runpt] << "\n" << 
                 m2vec[runpt] << "\n" << m3vec[runpt] << "\n" <<
                 m4vec[runpt] << "\n" << m5vec[runpt] << "\n" <<
                 m6vec[runpt] << "\n";
    }

    ndump.close();

    ndump.open("xvec.dat", std::ios_base::app);

    for (runpt=0; runpt<runpts; runpt++)
    {
        ndump << x_1vec[runpt] << "\n" << x0vec[runpt] << "\n" << x1vec[runpt] << "\n" << 
                 x2vec[runpt] << "\n" << x3vec[runpt] << "\n" <<
                 x4vec[runpt] << "\n" << x5vec[runpt] << "\n" <<
                 x6vec[runpt] << "\n";
    }

    ndump.close();

    ndump.open("cvec.dat", std::ios_base::app);

    for (runpt=0; runpt<runpts; runpt++)
    {
        ndump << c_1vec[runpt] << "\n" << c0vec[runpt] << "\n" << c1vec[runpt] << "\n" << 
                 c2vec[runpt] << "\n" << c3vec[runpt] << "\n" <<
                 c4vec[runpt] << "\n" << c5vec[runpt] << "\n" <<
                 c6vec[runpt] << "\n";
    }

    ndump.close();

    delete [] via;

}

int blocking::get_cell(int i, int j, int k, int dir_i, int dir_j, int dir_k, int b)
{
    int cell;

    if ((dir_i==1)&&(i==((L/b)-1)))
    {
        if ((dir_j==1)&&(j==((L/b)-1)))
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)(0*pow(L/b,2)+0*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(0*pow(L/b,2)+0*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(0*pow(L/b,2)+0*(L/b)+k);
            }
            else
            {
                cell = (int)(0*pow(L/b,2)+0*(L/b)+(k+dir_k));
            }
        }
        else if ((dir_j==(-1))&&(j==0))
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)(0*pow(L/b,2)+((L/b)-1)*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(0*pow(L/b,2)+((L/b)-1)*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(0*pow(L/b,2)+((L/b)-1)*(L/b)+k);
            }
            else
            {
                cell = (int)(0*pow(L/b,2)+((L/b)-1)*(L/b)+(k+dir_k));
            }
        }
        else if (dir_j==0)
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)(0*pow(L/b,2)+j*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(0*pow(L/b,2)+j*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(0*pow(L/b,2)+j*(L/b)+k);
            }
            else
            {
                cell = (int)(0*pow(L/b,2)+j*(L/b)+(k+dir_k));
            }
        }
        else
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)(0*pow(L/b,2)+(j+dir_j)*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(0*pow(L/b,2)+(j+dir_j)*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(0*pow(L/b,2)+(j+dir_j)*(L/b)+k);
            }
            else
            {
                cell = (int)(0*pow(L/b,2)+(j+dir_j)*(L/b)+(k+dir_k));
            }
        }
    }
    else if ((dir_i==(-1))&&(i==0))
    {
        if ((dir_j==1)&&(j=((L/b)-1)))
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+0*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+0*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+0*(L/b)+k);
            }
            else
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+0*(L/b)+(k+dir_k));
            }
        }
        else if ((dir_j==(-1))&&(j==0))
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+(L-1)*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+((L/b)-1)*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+((L/b)-1)*(L/b)+k);
            }
            else
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+((L/b)-1)*(L/b)+(k+dir_k));
            }
        }
        else if (dir_j==0)
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+j*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+j*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+j*(L/b)+k);
            }
            else
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+j*(L/b)+(k+dir_k));
            }
        }
        else
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+(j+dir_j)*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+(j+dir_j)*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+(j+dir_j)*(L/b)+k);
            }
            else
            {
                cell = (int)(((L/b)-1)*pow(L/b,2)+(j+dir_j)*(L/b)+(k+dir_k));
            }
        }
    }
    else if (dir_i==0)
    {
        if ((dir_j==1)&&(j=((L/b)-1)))
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)(i*pow(L/b,2)+0*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(i*pow(L/b,2)+0*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(i*pow(L/b,2)+0*(L/b)+k);
            }
            else
            {
                cell = (int)(i*pow(L/b,2)+0*(L/b)+(k+dir_k));
            }
        }
        else if ((dir_j==(-1))&&(j==0))
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)(i*pow(L/b,2)+((L/b)-1)*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(i*pow(L/b,2)+((L/b)-1)*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(i*pow(L/b,2)+((L/b)-1)*(L/b)+k);
            }
            else
            {
                cell = (int)(i*pow(L/b,2)+((L/b)-1)*(L/b)+(k+dir_k));
            }
        }
        else if (dir_j==0)
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)(i*pow(L/b,2)+j*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(i*pow(L/b,2)+j*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(i*pow(L/b,2)+j*(L/b)+k);
            }
            else
            {
                cell = (int)(i*pow(L/b,2)+j*(L/b)+(k+dir_k));
            }
        }
        else
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)(i*pow(L/b,2)+(j+dir_j)*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(i*pow(L/b,2)+(j+dir_j)*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(i*pow(L/b,2)+(j+dir_j)*(L/b)+k);
            }
            else
            {
                cell = (int)(i*pow(L/b,2)+(j+dir_j)*(L/b)+(k+dir_k));
            }
        }
    }
    else
    {
        if ((dir_j==1)&&(j=((L/b)-1)))
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+0*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+0*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+0*(L/b)+k);
            }
            else
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+0*(L/b)+(k+dir_k));
            }
        }
        else if ((dir_j==(-1))&&(j==0))
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+((L/b)-1)*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+((L/b)-1)*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+((L/b)-1)*(L/b)+k);
            }
            else
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+((L/b)-1)*(L/b)+(k+dir_k));
            }
        }
        else if (dir_j==0)
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+j*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+j*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+j*(L/b)+k);
            }
            else
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+j*(L/b)+(k+dir_k));
            }
        }
        else
        {
            if ((dir_k==1)&&(k==((L/b)-1)))
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+(j+dir_j)*(L/b)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+(j+dir_j)*(L/b)+((L/b)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+(j+dir_j)*(L/b)+k);
            }
            else
            {
                cell = (int)((i+dir_i)*pow(L/b,2)+(j+dir_j)*(L/b)+(k+dir_k));
            }
        }
    }

    return cell;
}

double blocking::separation(double r1[3], double r2[3])
{
    double d[3];
    d[0] = 0.0; d[1] = 0.0; d[2] = 0.0;

    double r;
    int dimR;

    for (dimR=0; dimR<dim; dimR++)
    {
        if ((r1[dimR]-r2[dimR])>(double)(L/2))
        {
            d[dimR] = r1[dimR]-r2[dimR] - L;
        }
        else if((r1[dimR]-r2[dimR])<(double)(-L/2))
        {
            d[dimR] = r1[dimR]-r2[dimR] + L;
        }
        else
        {
            d[dimR] = r1[dimR]-r2[dimR];
        }
    }

    r = sqrt(pow(d[0],2)+pow(d[1],2)+pow(d[2],2));

    return r;
}

blocking::~blocking()
{
}
