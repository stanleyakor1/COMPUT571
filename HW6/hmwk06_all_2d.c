
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

enum
{
    DIRICHLET=0,
    NEUMANN
};

enum
{
    LEFT=0,
    RIGHT,
    BOTTOM,
    TOP
};

int bc_type[4];

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

double u_true(double x,double y)
{
    return sin(2*x)*sin(2*y);
}

double ux_true(double x, double y)
{
    return 2*cos(2*x)*sin(2*y);
}

double uxx_true(double x, double y)
{
    return -4*u_true(x,y);
}

double uy_true(double x, double y)
{
    return 2*sin(2*x)*cos(2*y);
}

double uyy_true(double x, double y)
{
    return -4*u_true(x,y);
}

double un_true(int iside,double x,double y)
{
    double un;
    if (iside == LEFT)
        un = -ux_true(x,y);
    
    else if (iside == RIGHT)
        un = ux_true(x,y);
    
    else if (iside == BOTTOM)
        un = -uy_true(x,y);
    
    else if (iside == TOP)
        un = uy_true(x,y);

    return un;
}

double upp_true(double x,double y)
{
    return uxx_true(x,y) + uyy_true(x,y);
}
