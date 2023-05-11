
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* allocate_1d(int n)
{
    double *mem = (double*) malloc(n*sizeof(double));
    return &mem[0];
}

void free_1d(double **y)
{
    free(*y);
    *y = NULL;
}

double  f(double x)
{
    return cos(x);
}

int main(int argv, char** argc)
{
    int N = 64;
    
    double a = 0; 
    double b = 1; 
    
    double *x = allocate_1d(N+1);    
    double *y = allocate_1d(N+1);

    double h = (b-a)/N;
    
    for(int i = 0; i < N+1; i++)
    {
        x[i] = a + h*i;
        y[i] = f(x[i]);
    }
    
    double fxy[N+1][N+1];
    
    /* Assign valuels to f(x,y) using a double loop */
    
    
    /* Add statements to write out data x, y and metadata N, a and b.  */
    
    free_1d(&x);
    free_1d(&y);
}
