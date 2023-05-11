
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <math.h>

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

double sinc(double x,double y)
{
    double r = sqrt(pow(x,2) + pow(y,2));
    if (r == 0)
        return 1;
    else
        return sin(r)/r;
}

// # Test output using a shift.
double xshift = 3.0/8.0;
double yshift = 3.0/4.0;

double qtrue(double x, double y)
{
    double x1 = x - xshift;
    double y1 = y - yshift;
    //return (x1*x1 + y1*y1)/4.0;
    //return sin(2*M_PI*x1)*sin(2*M_PI*y1);
    return sinc(20*x1,20*y1);
}
