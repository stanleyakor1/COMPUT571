
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

double c2 = 0;
double c1 = -0.930763853398148;

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


double utrue(double x)
{
    double pi = M_PI;
    double utrue = (sqrt(pi)*x*erf(x) + exp(-x*x))/2.0 + c2*x + c1;
    return utrue;    
}

double rhs(double x)
{
    double upp = exp(-x*x);
    return upp;
}

void matvec(int N, double *u, double *L)
{
    // Apply BC
    u[-1] =  -u[1];
    
    u[N+1] = - u[N-1];
    
    for(int i = 0; i < N+1; i++)
    {
        L[i] = (u[i-1] - 2*u[i] + u[i+1]); 
    }
}

void jacobi(int N, double *F, double *u, double tol, int kmax, int prt)
{
    // Jacobi iteration
    double *u_old = allocate_1d(N+1, 1);
    double *Lu = allocate_1d(N+1,0);

    double update;
    int i, j, iter;
    double error;
  
    for (i = 0; i <= N; i++) {
        u_old[i] = u[i] = 0;
    }
    
    // Iterate using Jacobi method
    for (iter = 0; iter < kmax; iter++) 
    {
        matvec(N,u_old,Lu);
        error = 0.0;
        for (i = 0; i <= N; i++) {
            update= -0.5*(F[i] - Lu[i]);
            u[i] = u_old[i]+update ;
            error = fmax(fabs(update), error);
        }

        
       // printf("%5d %12.4e\n",iter,error);
        if (error < tol) 
        {
            break;
        }

        // Copy u to u_old
        for (i = 0; i <= N; i++) {
            u_old[i]=u[i];
        }
    }
    printf("Number of iteration = %d", iter);

    free_1d(&u_old,1);
    free_1d(&Lu,0);
}




int main(int argc, char** argv)
{
    int N = atoi(argv[1]);
    double *u = allocate_1d(N+1, 1);
    double *u_true = allocate_1d(N+1,1);
    double *F = allocate_1d(N+1,0);
    double *LU = allocate_1d(N+1, 0);
    int prt = 1;
    double a = -1.0;
    double b = 1.0;
    double h = (b - a)/N;
    double tol = 1e-12;
    int kmax = 1e6;
    for (int i =0; i<=N; i++)
    {
        F[i] = pow(h,2)*rhs(a + i*h);
        u_true[i] = utrue(a + i*h);
    } 
    u_true[-1] = -u_true[1];
    u_true[N+1] = -u_true[N-1];
    jacobi(N, F,u,tol, kmax, prt);
   
   matvec(N,u,LU);
    
   double er;
   double maxres = 0.0;
   double maxerr = 0.0;
    for(int i = 0; i < N+1; i++)
    {
        double r = fabs(F[i] - LU[i]);
        maxres = fmax(r,maxres);
    }
    for (int i = 0; i<=N; i++)
    {
        double er = fabs(u[i] - u_true[i]);
        maxerr = fmax(er,maxerr);
    }
    
    printf("\n");
    printf("Error norm : %12.4e\n",maxerr);
    printf("Residual : %12.4e\n",maxres);
    printf("-------------------------\n");
    
    FILE *fout = fopen("prob1.dat","w");
    if (fout == NULL) 
    {
        perror("Unable to create file");
        exit(1);
    }
    fwrite(&N, sizeof(int), 1, fout);
    fwrite(&a, sizeof(double), 1, fout);
    fwrite(&b, sizeof(double), 1, fout);
    fwrite(&u[0], sizeof(double), N+1, fout);
    fclose(fout);
    
    free_1d(&u_true,1);
    free_1d(&u,1);
    free_1d(&LU,0);
    free_1d(&F,0);
    return 0;
    
}

