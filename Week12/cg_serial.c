
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

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

double utrue(double x)
{
    return x*(x-1)/2;
    //return sin(x);    
}

double rhs(double x)
{
    return 1;
    //return -sin(x);
}

void matvec(int N, double *u, double *L)
{
    for(int i = 0; i < N+1; i++)
        L[i] = (u[i-1] - 2*u[i] + u[i+1]); 
}

int cg(int N, double *F, double *u, double tol, int kmax, int prt)
{
    double *uk = allocate_1d(N+1,1);
    double *pk = allocate_1d(N+1,1);
    
    double rk[N+1];
    for(int i = 0; i < N+1; i++)
    {
        uk[i] = 0;
        rk[i] = F[i];
        pk[i] = rk[i];    // # Start with uk = 0 --> r = b - Au = b
    }
                
    double Apk[N+1];
    double rkp1[N+1];
    
    int itcount = 0;
    for(int k = 0; k < kmax; k++)
    {
        matvec(N,pk,Apk);
        
        // # Compute the residual
        double rTr = 0;
        double pTAp = 0;
        for(int i = 0; i < N+1; i++)
        {
            rTr += rk[i]*rk[i];
            pTAp += pk[i]*Apk[i];
        }
        double alpha = rTr/pTAp;
        
        double rpTrp = 0;
        for(int i = 0; i < N+1; i++)
        {
            uk[i] = uk[i] + alpha*pk[i];
            rkp1[i] = rk[i] - alpha*Apk[i];
            rpTrp += rkp1[i]*rkp1[i];
        }
        
        double beta = rpTrp/rTr;
        
        // update search direction
        for(int i = 0; i < N+1; i++)
            pk[i] = rkp1[i] + beta*pk[i];
        
        if (prt)
            printf("%5d %12.4e\n",k,rpTrp);
        
        itcount = k;
        if (rpTrp < tol)
            break;
        
        // # Update values
        for(int i = 0; i < N+1; i++)
            rk[i] = rkp1[i];
        
        rTr = rpTrp;
            
    }
    for(int i = 0; i < N+1; i++)
        u[i] = uk[i];    
    
    return itcount;
}

int main(int argc, char** argv)
{

    int N;
    if (argc < 2)
        N = 64;
    else
        N = atoi(argv[1]);
    
    
    // # Physical parameters
    double a = 0; 
    double b = 1;
    
    // # Numerial parameters
    double tol = 1e-12;
    int kmax = 100000;
    int prt = 0;
    double h = (b-a)/N;
    
    // # Arrays
    double *u = allocate_1d(N+1,1);
    double F[N+1];
    double Lu[N+1];
    
    // # Initialization
    for(int i = 0; i < N+1; i++)
    {
        double x =  a + i*h;
        F[i] = h*h*rhs(x);        
        u[i] = 0;
    }
    
    // # Compute the solution using Jacobi method
    int itcount = cg(N,F,u,tol,kmax,prt);
        
    
    /* Supply correct ghost cells values */
    matvec(N,u,Lu);
    
    // # Compute error and residual
    double max_res = 0;
    
    for(int i = 0; i < N+1; i++)
    {
        double x = a + h*i;
        
        double res = fabs(F[i] - Lu[i]);
        max_res = fmax(res,max_res);            
        
#if 0
        if (N < 32)
            printf("%5d %16.8f %12.4e %12.4e\n",i, u[i],err,res);
#endif        
    }
    
    printf("\n");
    printf("N = %d\n",N);
    printf("Iteration count     : %24d\n",itcount);
    printf("Residual (inf-norm) : %24.4e\n",max_res/2);
    
    FILE *fout = fopen("cg.dat","wb");   
    fwrite(&N,sizeof(int),1,fout);
    fwrite(&a,sizeof(double),1,fout);
    fwrite(&b,sizeof(double),1,fout);
    fwrite(&u[0],sizeof(double),N+1, fout); 
    fclose(fout);
            
    return 0;
}    
