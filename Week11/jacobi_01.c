
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

double* allocate_1d(int n, int m)
{
    double *mem = (double*) malloc((n + 2*m)*sizeof(double));
    return &mem[m];
}

void free_1d(double **x, int m)
{
    free(&(*x)[-m]);
    *x = NULL;
}


void get_rhs(int N, double a, double b, double *F)
{
    double h = (b-a)/N;
    for(int i = 0; i < N+1; i++)
    {
        F[i] = h*h;        
    }    
}

double utrue(double x)
{
    double utrue = x*(x-1)/2;
    
    return utrue;
}

void apply_laplacian(int N, double *u, double *L)
{
    double g0 = 0;
    u[-1] = 2*g0 -u[1];
    
    double g1 = 0;
    u[N+1] = 2*g1 - u[N-1];
    
    for(int i = 0; i<N+1; i++)
        L[i] = (u[i-1] - 2*u[i] + u[i+1]);
}

void jacobi_ver2(int N, double *F, double *u, double tol, int kmax, int prt)
{
    for (int j = 0; j<N+1; j++)
        u[j] = 0;
    /* Jacobi iteration */
    double *Lu = allocate_1d(N+1,0);
    double *uk = allocate_1d(N+1,0);
    double *ukp1 = allocate_1d(N+1,1);
    for (int k = 0; k<kmax; k++)
    {

        apply_laplacian(N, a, b, u, Lu);
        double di = -2;
#if 0
        
#endif
    }
}

int main(int argc, char** argv)
{
    
    double a = 0; 
    double b = 1;
    
    int N = 8;
    double tol = 1e-12;
    int kmax = 1000;
    int prt = 1;
        
    double *u = allocate_1d(N,1);
    double *F = allocate_1d(N,0);
    
    get_rhs(N,a,b,F);
    
    jacobi_ver2(N,F,u,tol,kmax,prt);
        
    // # TODO : Compute the error
    
    // # TODO : Compute the residual
    
    
    return 0;
}    
