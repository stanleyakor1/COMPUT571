
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hmwk06_all_1d.c"

int bc_type[2];
void matvec(int N, double *u, double *L)
{
    // # TODO : Generalize boundary condtions. 
    int i1 =0;
    int i2 = N;
     
    u[-1] = ( bc_type[LEFT] == NEUMANN)? u[0]:-u[0];
    u[N] = (bc_type[RIGHT] == NEUMANN)? u[N-1]: -u[N-1];
    
    for(int i = i1; i < i2; i++)
    {
      L[i] = u[i-1] - 2*u[i] + u[i+1];
    }
  
}

int cg(int N, double *F, double *u, double tol, int kmax, int prt)
{
    
    // # TODO : Check start and stop indices. 
    
    int i1 = 0;
    int i2 = N;
    
    double *uk = allocate_1d(N,0);
    double *pk = allocate_1d(N,1);    
    double *rk = allocate_1d(N,0);
    double *Apk = allocate_1d(N,0);
    double *rkp1 = allocate_1d(N,0);
    
    
    for(int i = i1; i < i2; i++)
    {
        uk[i] = 0;
        rk[i] = F[i];
        pk[i] = rk[i];    // # Start with uk = 0 --> r = b - Au = b            
    }
    
    int itcount = 0;
    for(int k = 0; k < kmax; k++)
    {        
        matvec(N,pk,Apk);
        
        double rTr = 0;
        double pTAp = 0;
        for(int i = i1; i < i2; i++)
        {
            rTr += rk[i]*rk[i];
            pTAp += pk[i]*Apk[i];
        }
        if (pTAp == 0)
        {
            printf("pTAp == 0; returning solution\n");
            break;
        }
        double alpha = rTr/pTAp;
        
        double rpTrp = 0;
        double max_res = 0;
        for(int i = i1; i < i2; i++)
        {
            uk[i] = uk[i] + alpha*pk[i];
            rkp1[i] = rk[i] - alpha*Apk[i];
            rpTrp += rkp1[i]*rkp1[i];
            max_res = fmax(fabs(rkp1[i]),max_res);
        }

        itcount = k+1;
        if (prt)
            printf("%5d %12.4e\n",itcount,max_res);
        
        if (max_res < tol)
            break;
                
        double beta = rpTrp/rTr;
        
        // # Update search directions
        for(int i = i1; i < i2; i++)
            pk[i] = rkp1[i] + beta*pk[i];
                
        // # Update residuals
        for(int i = i1; i < i2; i++)
            rk[i] = rkp1[i];
        
        rTr = rpTrp;
            
    }
    for(int i = i1; i < i2; i++)
        u[i] = uk[i];    
    
    free_1d(&uk,0);
    free_1d(&pk,1);
    free_1d(&rk,0);
    free_1d(&Apk,0);
    free_1d(&rkp1,0);
        
    return itcount;
}

int main(int argc, char** argv)
{

    int N;
    int code = parse_input(argc,argv,&N,&bc_type[0]);
    if (code > 0)
        exit(1);
    
    // # Domain
    double a = 0, b = 1;
    
    // # Numerical parameters
    int kmax;
    if (bc_type[LEFT] == NEUMANN && bc_type[RIGHT] == NEUMANN)
    {
        kmax = 30;
        printf("All Neumann boundary conditions may not converge. " \
               "Value of kmax is set to %d.\n\n",kmax);
    }
    else
        kmax = 100000;
        
    double tol = 1e-12;
    int prt = 0;

    double h = (b-a)/N;
    
    // # Arrays    
    double *u = allocate_1d(N,1);
    double *F = allocate_1d(N,0);
    double *x = allocate_1d(N,0);
    
    // # Initialization
    for(int i = 0; i < N; i++)
    {
        x[i] = a + (i+0.5)*h;
        F[i] = h*h*upp_true(x[i]);        
    }

    // # TODO : Modify right hand side to account for inhomogeneous BCs
    
    if (bc_type[LEFT] == DIRICHLET)
        F[0] -=2*u_true(a);
    
    if (bc_type[RIGHT] == DIRICHLET)
        F[N-1] -=2*u_true(b);
    
    if (bc_type[LEFT] == NEUMANN)
        F[0] -=h*up_true(LEFT,a);
    
    if (bc_type[RIGHT] == NEUMANN)
        F[N-1] -=h*up_true(RIGHT,b);
    
    
    // # Compute the solution using conjugate Gradient method
    int itcount = cg(N,F,u,tol,kmax,prt);
            
    // # Compute error 
    double max_err = 0;
    for(int i = 0; i < N; i++)
        max_err = fmax(fabs(u_true(x[i]) - u[i]),max_err);
    
    // # Print out results
    printf("\n");
    printf("N = %d\n",N);
    printf("Iteration count : %12d\n",itcount);
    printf("Error (pde)     : %12.4e\n",max_err);        
    
    FILE* fout = fopen("err_02.dat","wb");   
    fwrite(&max_err,sizeof(double),1,fout);
    fwrite(&itcount,sizeof(int),1,fout);
    fclose(fout);    
    
    fout = fopen("cg_02.dat","wb");   
    fwrite(&N,sizeof(int),1,fout);
    fwrite(bc_type,sizeof(int),2,fout);
    fwrite(&a,sizeof(double),1,fout);
    fwrite(&b,sizeof(double),1,fout);
    fwrite(&x[0],sizeof(double),N, fout); 
    fwrite(&u[0],sizeof(double),N, fout);             
    fclose(fout);

    free_1d(&u,1);
    free_1d(&F,0);
    free_1d(&x,0);
    
    return 0;
}    
