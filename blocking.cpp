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
    lattice.set_dimensions(1);
    L   = lattice.L; // length of lattice (number of sites)
    dim = lattice.dim; // dimensionality of lattice
    nL  = pow(L,dim);
}

blocking::blocking(int b)
{
    coarse_lattice.set_dimensions(b);
    Lb = coarse_lattice.L;

    //lattice_map(L, b, bshift);
    //lattice_map.U_L_Lb; // the coarse cell that each original lattice site belongs to

    int n_configs;
    n_configs = 20;

    std::ostringstream filenameStream;
    std:string filename;

    ifstream fin;
    string line;

    std::ofstream ndump;

    int i,j,k;
    int ib,jb,kb;

    // read in data

    int n[nL];

    // coarse-grained n (average per cell)
    double nb[(int)pow(Lb,dim)];

    // collective variables for each MC config
    double n_1; // sum\, n_i
    double n0;  // sum_NN\, n_i * n_i+1 (i=x,y,z)
    double n1;  // sum_{diagonal in plane}\, n_i * n_j
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

    double n_1vec[n_configs];
    double n0vec[n_configs];
    double n1vec[n_configs];
    double n2vec[n_configs];
    double n3vec[n_configs];
    double n4vec[n_configs];
    double n5vec[n_configs];
    double n6vec[n_configs];

    int celli,cellj,cellk,celll;

    double nn;

    double nn_i,nn_j,nn_k,nn_l;

    for (i=0; i<(int)(pow(Lb,dim)); i++)
    {
        nb[i] = 0.0;
    }

    int config;

    for (config=0; config<n_configs; config++)
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

        filenameStream << "./netlib/n" << config << ".dat";
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

        i = 0; j = 0; k = 0;

        ib = 0; jb = 0; kb = 0;

        for (i=0; i<(int)(pow(Lb,dim)); i++)
        {
            nb[i] = 0.0;
        }

        n_1 = 0.0; n0 = 0.0; n1 = 0.0; n2 = 0.0; n3 = 0.0; n4 = 0.0; n5 = 0.0; n6 = 0.0;

        nn = 0.0;

        nn_i = 0.0; nn_j = 0.0; nn_k = 0.0; nn_l = 0.0;

        celli = 0; cellj = 0; cellk = 0; celll = 0;

    #pragma omp parallel for shared (n,nb) private (i,j,k,celli,nn,ib,jb,kb)
        for (i=0; i<Lb; i++)
        {
            for (j=0; j<Lb; j++)
            {
                for (k=0; k<Lb; k++)
                {
                    // calculate average n for this cell block

                    celli = (int)(i*pow(Lb,2)+j*(Lb)+k);

        #pragma omp atomic write
                    nb[celli] = 0.0;

                    nn = 0.0;

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
                                        }
                                        else
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb-L)*L+k*b+kb)];
                                        }
                                    }
                                    else
                                    {
                                        if ((k*b+kb)>(L-1))
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb)*L+k*b+kb-L)];
                                        }
                                        else
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib-L)*pow(L,2)+(j*b+jb)*L+k*b+kb)];
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
                                        }
                                        else
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib)*pow(L,2)+(j*b+jb-L)*L+k*b+kb)];

                                        }
                                    }
                                    else
                                    {
                                        if ((k*b+kb)>(L-1))
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib)*pow(L,2)+(j*b+jb)*L+k*b+kb-L)];
                                        }
                                        else
                                        {
        #pragma omp atomic read
                                            nn = n[(int)((i*b+ib)*pow(L,2)+(j*b+jb)*L+k*b+kb)];
                                        }
                                    }
                                }

        #pragma omp atomic write
                                nb[celli] = nb[celli] + nn;

                            }
                        }
                    }
                }
            }
        }

    #pragma omp parallel for shared (n,nb) private (i,j,k,celli,nn) reduction(+:n_1)
        for (i=0; i<Lb; i++)
        {
            for (j=0; j<Lb; j++)
            {
                for (k=0; k<Lb; k++)
                {
                    celli = (int)(i*pow(Lb,2)+j*(Lb)+k);

        #pragma omp atomic read
                    nn = nb[celli];

        #pragma omp atomic write
                    nb[celli] = nn/(double)(pow(b,dim));

                    n_1 += nn/(double)(pow(b,dim));
                }
            }
        }

        ndump.open("nb.dat", std::ios_base::app);

        for (i=0; i<(int)(pow(Lb,dim)); i++)
        {
            ndump << nb[i] << "\n";
        }
        ndump.close();

    #pragma omp parallel for shared (nb) private (i,j,k,celli,cellj,nn_i,nn_j) reduction(+:n0)
        for (i=0; i<Lb; i++)
        {
            for (j=0; j<Lb; j++)
            {
                for (k=0; k<Lb; k++)
                {
                    // sum_NN\, n_i * n_i+1 (i=x,y,z)

                    celli = (int)(i*pow(Lb,2)+j*(Lb)+k);
        #pragma omp atomic read
                    nn_i = nb[celli];

                    cellj = get_cell(i,j,k,1,0,0,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n0 += nn_i*nn_j;

                    cellj = get_cell(i,j,k,0,1,0,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n0 += nn_i*nn_j;

                    cellj = get_cell(i,j,k,0,0,1,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n0 += nn_i*nn_j;
                }
            }
        }

    #pragma omp parallel for shared (nb) private (i,j,k,celli,cellj,nn_i,nn_j) reduction(+:n1)
        for (i=0; i<Lb; i++)
        {
            for (j=0; j<Lb; j++)
            {
                for (k=0; k<Lb; k++)
                {
                    // sum_{diagonal in plane}\, n_i * n_j
                    // e.g. n_(0,0,0) * n_(1,1,0)

                    celli = (int)(i*pow(Lb,2)+j*(Lb)+k);
        #pragma omp atomic read
                    nn_i = nb[celli];

                    cellj = get_cell(i,j,k,1,1,0,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n1 += nn_i*nn_j; 

                    cellj = get_cell(i,j,k,1,-1,0,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n1 += nn_i*nn_j; 

                    cellj = get_cell(i,j,k,1,0,1,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n1 += nn_i*nn_j; 

                    cellj = get_cell(i,j,k,1,0,-1,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n1 += nn_i*nn_j; 

                    cellj = get_cell(i,j,k,0,1,1,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n1 += nn_i*nn_j; 

                    cellj = get_cell(i,j,k,0,1,-1,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n1 += nn_i*nn_j; 
                }
            }
        }

    #pragma omp parallel for shared (nb) private (i,j,k,celli,cellj,nn_i,nn_j) reduction(+:n2)
        for (i=0; i<Lb; i++)
        {
            for (j=0; j<Lb; j++)
            {
                for (k=0; k<Lb; k++)
                {
                    // sum_{cubic diagonal}\, n_i * n_j
                    // e.g. n_(0,0,0) * n_(1,1,1)

                    celli = (int)(i*pow(Lb,2)+j*(Lb)+k);
        #pragma omp atomic read
                    nn_i = nb[celli];

                    cellj = get_cell(i,j,k,1,1,1,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n2 += nn_i*nn_j; 

                    cellj = get_cell(i,j,k,1,1,-1,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n2 += nn_i*nn_j; 

                    cellj = get_cell(i,j,k,1,-1,-1,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n2 += nn_i*nn_j; 

                    cellj = get_cell(i,j,k,1,-1,1,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n2 += nn_i*nn_j; 
                }
            }
        }

    #pragma omp parallel for shared (nb) private (i,j,k,celli,cellj,cellk,celll,nn,nn_i,nn_j,nn_k,nn_l) reduction(+:n3)
        for (i=0; i<Lb; i++)
        {
            for (j=0; j<Lb; j++)
            {
                for (k=0; k<Lb; k++)
                {
                    // sum_{principal planes}\, n_i * n_j * n_k * n_l
                    // e.g. n_(0,0,0) * n_(0,1,0) * n_(1,0,0) * n_(1,1,0)

                    celli = (int)(i*pow(Lb,2)+j*(Lb)+k);
        #pragma omp atomic read
                    nn_i = nb[celli];

                    cellj = get_cell(i,j,k,1,0,0,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    cellk = get_cell(i,j,k,0,1,0,b);

        #pragma omp atomic read
                    nn_k = nb[cellk];

                    celll = get_cell(i,j,k,1,1,0,b);

        #pragma omp atomic read
                    nn_l = nb[celll];

                    n3 += nn_i*nn_j*nn_k*nn_l; 

                    cellj = get_cell(i,j,k,1,0,0,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    cellk = get_cell(i,j,k,0,0,1,b);

        #pragma omp atomic read
                    nn_k = nb[cellk];

                    celll = get_cell(i,j,k,1,0,1,b);

        #pragma omp atomic read
                    nn_l = nb[celll];

                    n3 += nn_i*nn_j*nn_k*nn_l; 

                    cellj = get_cell(i,j,k,0,1,0,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    cellk = get_cell(i,j,k,0,0,1,b);

        #pragma omp atomic read
                    nn_k = nb[cellk];

                    celll = get_cell(i,j,k,0,1,1,b);

        #pragma omp atomic read
                    nn_l = nb[celll];

                    n3 += nn_i*nn_j*nn_k*nn_l; 
                }
            }
        }

    #pragma omp parallel for shared (nb) private (i,j,k,celli,cellj,cellk,celll,nn,nn_i,nn_j,nn_k,nn_l) reduction(+:n4)
        for (i=0; i<Lb; i++)
        {
            for (j=0; j<Lb; j++)
            {
                for (k=0; k<Lb; k++)
                {
                    // sum_{diagonal planes}\, n_i * n_j * n_k * n_l
                    // e.g. n_(0,0,0) * n_(1,0,0) * n_(0,1,1) * n_(1,1,1)

                    celli = (int)(i*pow(Lb,2)+j*(Lb)+k);
        #pragma omp atomic read
                    nn_i = nb[celli];

                    cellj = get_cell(i,j,k,1,0,0,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    cellk = get_cell(i,j,k,0,1,1,b);

        #pragma omp atomic read
                    nn_k = nb[cellk];

                    celll = get_cell(i,j,k,1,1,1,b);

        #pragma omp atomic read
                    nn_l = nb[celll];

                    n4 += nn_i*nn_j*nn_k*nn_l; 

                    cellj = get_cell(i,j,k,1,0,0,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    cellk = get_cell(i,j,k,0,1,-1,b);

        #pragma omp atomic read
                    nn_k = nb[cellk];

                    celll = get_cell(i,j,k,1,1,-1,b);

        #pragma omp atomic read
                    nn_l = nb[celll];

                    n4 += nn_i*nn_j*nn_k*nn_l; 

                    cellj = get_cell(i,j,k,0,1,0,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    cellk = get_cell(i,j,k,1,0,1,b);

        #pragma omp atomic read
                    nn_k = nb[cellk];

                    celll = get_cell(i,j,k,1,1,1,b);

        #pragma omp atomic read
                    nn_l = nb[celll];

                    n4 += nn_i*nn_j*nn_k*nn_l; 

                    cellj = get_cell(i,j,k,0,1,0,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    cellk = get_cell(i,j,k,1,0,-1,b);

        #pragma omp atomic read
                    nn_k = nb[cellk];

                    celll = get_cell(i,j,k,1,1,-1,b);

        #pragma omp atomic read
                    nn_l = nb[celll];

                    n4 += nn_i*nn_j*nn_k*nn_l; 

                    cellj = get_cell(i,j,k,0,0,1,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    cellk = get_cell(i,j,k,1,1,0,b);

        #pragma omp atomic read
                    nn_k = nb[cellk];

                    celll = get_cell(i,j,k,1,1,1,b);

        #pragma omp atomic read
                    nn_l = nb[celll];

                    n4 += nn_i*nn_j*nn_k*nn_l; 

                    cellj = get_cell(i,j,k,0,0,1,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    cellk = get_cell(i,j,k,1,-1,0,b);

        #pragma omp atomic read
                    nn_k = nb[cellk];

                    celll = get_cell(i,j,k,1,-1,1,b);

        #pragma omp atomic read
                    nn_l = nb[celll];

                    n4 += nn_i*nn_j*nn_k*nn_l; 
                }
            }
        }

    #pragma omp parallel for shared (nb) private (i,j,k,celli,cellj,cellk,celll,nn,nn_i,nn_j,nn_k,nn_l) reduction(+:n5)
        for (i=0; i<Lb; i++)
        {
            for (j=0; j<Lb; j++)
            {
                for (k=0; k<Lb; k++)
                {
                    // sum_{tetrahedral vertices}\, n_i * n_j * n_k * n_l
                    // e.g. n_(0,0,0) * n_(1,0,1) * n_(0,1,1) * n_(1,1,0)

                    celli = (int)(i*pow(Lb,2)+j*(Lb)+k);
        #pragma omp atomic read
                    nn_i = nb[celli];

                    cellj = get_cell(i,j,k,1,0,1,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    cellk = get_cell(i,j,k,0,1,1,b);

        #pragma omp atomic read
                    nn_k = nb[cellk];

                    celll = get_cell(i,j,k,1,1,0,b);

        #pragma omp atomic read
                    nn_l = nb[celll];

                    n5 += nn_i*nn_j*nn_k*nn_l; 

                    cellj = get_cell(i,j,k,1,0,1,b);

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    cellk = get_cell(i,j,k,0,-1,1,b);

        #pragma omp atomic read
                    nn_k = nb[cellk];

                    celll = get_cell(i,j,k,1,-1,0,b);

        #pragma omp atomic read
                    nn_l = nb[celll];

                    n5 += nn_i*nn_j*nn_k*nn_l; 
                }
            }
        }

    #pragma omp parallel for shared (nb) private (i,j,k,celli,cellj,nn_i,nn_j) reduction(+:n6)
        for (i=0; i<Lb; i++)
        {
            for (j=0; j<Lb; j++)
            {
                for (k=0; k<Lb; k++)
                {
                    // sum_NNN\, n_i * n_i+2 (i=x,y,z)

                    celli = (int)(i*pow(Lb,2)+j*(Lb)+k);
        #pragma omp atomic read
                    nn_i = nb[celli];

                    //cellj = get_cell(i,j,k,2,0,0,b); 
                    if (i==((Lb)-1))
                    {
                        cellj = get_cell(0,j,k,1,0,0,b);
                    }
                    else
                    {
                        cellj = get_cell(i+1,j,k,1,0,0,b);
                    }

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n6 += nn_i*nn_j; 

                    //cellj = get_cell(i,j,k,0,2,0,b); 
                    if (j==((Lb)-1))
                    {
                        cellj = get_cell(i,0,k,0,1,0,b);
                    }
                    else
                    {
                        cellj = get_cell(i,j+1,k,0,1,0,b);
                    }

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n6 += nn_i*nn_j; 

                    //cellj = get_cell(i,j,k,0,0,2,b); 
                    if (k==((Lb)-1))
                    {
                        cellj = get_cell(i,j,0,0,0,1,b);
                    }
                    else
                    {
                        cellj = get_cell(i,j,k+1,0,0,1,b);
                    }

        #pragma omp atomic read
                    nn_j = nb[cellj];

                    n6 += nn_i*nn_j; 
                }
            }
        }

        n_1vec[config] = n_1;
        n0vec[config] = n0;
        n1vec[config] = n1;
        n2vec[config] = n2;
        n3vec[config] = n3;
        n4vec[config] = n4;
        n5vec[config] = n5;
        n6vec[config] = n6;

    }

    // save results to disk

    ndump.open("nvec.dat", std::ios_base::app);

    for (config=0; config<n_configs; config++)
    {
        ndump << n_1vec[config] << "\n" << n0vec[config] << "\n" << n1vec[config] << "\n" << 
                 n2vec[config] << "\n" << n3vec[config] << "\n" <<
                 n4vec[config] << "\n" << n5vec[config] << "\n" <<
                 n6vec[config] << "\n";
    }

    ndump.close();

}

int blocking::get_cell(int i, int j, int k, int dir_i, int dir_j, int dir_k, int b)
{
    int cell;

    if ((dir_i==1)&&(i==((Lb)-1)))
    {
        if ((dir_j==1)&&(j==((Lb)-1)))
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)(0*pow(Lb,2)+0*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(0*pow(Lb,2)+0*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(0*pow(Lb,2)+0*(Lb)+k);
            }
            else
            {
                cell = (int)(0*pow(Lb,2)+0*(Lb)+(k+dir_k));
            }
        }
        else if ((dir_j==(-1))&&(j==0))
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)(0*pow(Lb,2)+((Lb)-1)*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(0*pow(Lb,2)+((Lb)-1)*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(0*pow(Lb,2)+((Lb)-1)*(Lb)+k);
            }
            else
            {
                cell = (int)(0*pow(Lb,2)+((Lb)-1)*(Lb)+(k+dir_k));
            }
        }
        else if (dir_j==0)
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)(0*pow(Lb,2)+j*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(0*pow(Lb,2)+j*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(0*pow(Lb,2)+j*(Lb)+k);
            }
            else
            {
                cell = (int)(0*pow(Lb,2)+j*(Lb)+(k+dir_k));
            }
        }
        else
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)(0*pow(Lb,2)+(j+dir_j)*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(0*pow(Lb,2)+(j+dir_j)*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(0*pow(Lb,2)+(j+dir_j)*(Lb)+k);
            }
            else
            {
                cell = (int)(0*pow(Lb,2)+(j+dir_j)*(Lb)+(k+dir_k));
            }
        }
    }
    else if ((dir_i==(-1))&&(i==0))
    {
        if ((dir_j==1)&&(j=((Lb)-1)))
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+0*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+0*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+0*(Lb)+k);
            }
            else
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+0*(Lb)+(k+dir_k));
            }
        }
        else if ((dir_j==(-1))&&(j==0))
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+(L-1)*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+((Lb)-1)*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+((Lb)-1)*(Lb)+k);
            }
            else
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+((Lb)-1)*(Lb)+(k+dir_k));
            }
        }
        else if (dir_j==0)
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+j*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+j*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+j*(Lb)+k);
            }
            else
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+j*(Lb)+(k+dir_k));
            }
        }
        else
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+(j+dir_j)*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+(j+dir_j)*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+(j+dir_j)*(Lb)+k);
            }
            else
            {
                cell = (int)(((Lb)-1)*pow(Lb,2)+(j+dir_j)*(Lb)+(k+dir_k));
            }
        }
    }
    else if (dir_i==0)
    {
        if ((dir_j==1)&&(j=((Lb)-1)))
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)(i*pow(Lb,2)+0*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(i*pow(Lb,2)+0*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(i*pow(Lb,2)+0*(Lb)+k);
            }
            else
            {
                cell = (int)(i*pow(Lb,2)+0*(Lb)+(k+dir_k));
            }
        }
        else if ((dir_j==(-1))&&(j==0))
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)(i*pow(Lb,2)+((Lb)-1)*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(i*pow(Lb,2)+((Lb)-1)*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(i*pow(Lb,2)+((Lb)-1)*(Lb)+k);
            }
            else
            {
                cell = (int)(i*pow(Lb,2)+((Lb)-1)*(Lb)+(k+dir_k));
            }
        }
        else if (dir_j==0)
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)(i*pow(Lb,2)+j*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(i*pow(Lb,2)+j*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(i*pow(Lb,2)+j*(Lb)+k);
            }
            else
            {
                cell = (int)(i*pow(Lb,2)+j*(Lb)+(k+dir_k));
            }
        }
        else
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)(i*pow(Lb,2)+(j+dir_j)*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)(i*pow(Lb,2)+(j+dir_j)*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)(i*pow(Lb,2)+(j+dir_j)*(Lb)+k);
            }
            else
            {
                cell = (int)(i*pow(Lb,2)+(j+dir_j)*(Lb)+(k+dir_k));
            }
        }
    }
    else
    {
        if ((dir_j==1)&&(j=((Lb)-1)))
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+0*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+0*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+0*(Lb)+k);
            }
            else
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+0*(Lb)+(k+dir_k));
            }
        }
        else if ((dir_j==(-1))&&(j==0))
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+((Lb)-1)*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+((Lb)-1)*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+((Lb)-1)*(Lb)+k);
            }
            else
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+((Lb)-1)*(Lb)+(k+dir_k));
            }
        }
        else if (dir_j==0)
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+j*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+j*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+j*(Lb)+k);
            }
            else
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+j*(Lb)+(k+dir_k));
            }
        }
        else
        {
            if ((dir_k==1)&&(k==((Lb)-1)))
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+(j+dir_j)*(Lb)+0);
            }
            else if ((dir_k==(-1))&&(k==0))
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+(j+dir_j)*(Lb)+((Lb)-1));
            }
            else if (dir_k==0)
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+(j+dir_j)*(Lb)+k);
            }
            else
            {
                cell = (int)((i+dir_i)*pow(Lb,2)+(j+dir_j)*(Lb)+(k+dir_k));
            }
        }
    }

    return cell;
}


blocking::~blocking()
{
}
