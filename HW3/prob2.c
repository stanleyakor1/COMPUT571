
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

double f(double x)
{
    return cos(x)*exp(-10*x*x);
}

double fp(double x)
{
    /* Put true derivative here */
    return -exp(-10*x*x)*(sin(x)+20*x*cos(x));
}



double* allocate_1d(int n, int m)
{
    /* Use malloc to allocate memory */
    double *y = (double*) malloc((n+2*m)*sizeof(double));
    return y;
}


void free_1d(double **x, int m)
{
    /* Use free to free memory;  Set value of pointer to NULL after freeing memory.  */
    if (*x != NULL) 
    {
        free(*x+m);
        *x = NULL;
    }
    
}

int main(int argv, char**argc)
{
    int N = 8;
    int m = 1;
    double a = -1;
    double b = 1;
    double g;
    
    double *x = allocate_1d(N+1,m); 
    double h = (b-a)/N;

    for(int i = -m; i < N+m+1; i++)
    {
        x[i] = a + i*h;
    }
    printf("%6s %8s %12s %12s %12s\n","i","x","f'(x)","g(x)","err(x)");
    printf("----------------------------------------------------------------\n");
    for(int i = 0; i <= N+m+1; i++)
    {
        g = (f(x[i+1])-f(x[i-1]))/(2*h);
        if (i >= 0 && i < N+m)
            printf("%5d %12.8f %12.8f %12.8f %12.4e\n",i,x[i],fp(x[i]),g,fabs(g-fp(x[i])));
        
    }
    
    free_1d(&x,0);
    return 0;
}
