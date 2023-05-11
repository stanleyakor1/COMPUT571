
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
double initial_condition(double x, double y)
{
    double y0 = 0.5;
    double r0 = 0.2;
    double r = fabs(y - y0);    
    if (r <= r0)
        return 0.25*(1 + cos(M_PI*r/r0));
    else
        return 0;
}

void velocity(double x, double y, double *u, double *v)
{
    *u = pow(sin(M_PI*x),2)*sin(2*M_PI*y);
    *v = -pow(sin(M_PI*y),2)*sin(2*M_PI*x);
}
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

enum
{
    DIR_X = 0,
    DIR_Y
};
