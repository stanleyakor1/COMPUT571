
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <math.h>

enum
{
    LEFT=0,
    RIGHT,
    DOWN,
    UP
};

enum
{
    DIR_X = 0,
    DIR_Y
};

double* allocate_1d(int n, int m)
{
    double *mem = (double*) malloc((n + 2*m)*sizeof(double));
    return mem+m;
}

void free_1d(double **x, int m)
{
    free(*x-m);
    *x = NULL;
}

double** allocate_2d(int n, int m, int mbc)
{
    int rows = n + 2*mbc;
    int cols = m + 2*mbc; 

    double   *qmem = malloc(rows*cols*sizeof(double));
    double **qrows = malloc(rows*sizeof(double*));

    for(int i = 0; i < rows; i++)
    {
        qrows[i] = &qmem[cols*i + mbc];
    }    
    return &qrows[mbc];
}

void free_2d(double ***q,int mbc)
{
    free(&(*q)[-mbc][-mbc]);
    free(&(*q)[-mbc]);
    *q = NULL;
}

double qtrue(double x, double y)
{
    return (x*x + y*y)/4.0;
    //return sin(x)*sin(y);
}
