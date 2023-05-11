
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hmwk06_all_2d.c"



void matvec(int N, double **u, double **L)
{
    
    
    for(int j = 0; j < N; j++)
    {
        if (bc_type[LEFT] == DIRICHLET)
            u[-1][j] = -u[0][j];
        else
            u[-1][j] = u[0][j];
        if (bc_type[RIGHT] == DIRICHLET)
            u[N][j] = -u[N-1][j];
        else
            u[N][j] = u[N-1][j];
    }
    
    for(int j = 0; j < N; j++)
    {
        if (bc_type[BOTTOM] == DIRICHLET)
            u[j][-1] = -u[j][0];
        else
            u[j][-1] = u[j][0];
        
        if (bc_type[TOP] == DIRICHLET)
            u[j][N] = -u[j][N-1];
        else
            u[j][N] = u[j][N-1];
    }
    
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j <N; j++)
        {
            L[i][j] = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]- 4*u[i][j]) ;
        }
    }
    
    
}



int cg(int N, double **F, double **u, double tol, int kmax, int prt)
{
    // # TODO : Implement at 2d CG solver

    int itcount = 0;
    int i1 =0;
    int i2 = N;
    
    double **uk = allocate_2d(N,N,0);
    double **pk = allocate_2d(N,N,1);    
    double **rk = allocate_2d(N,N,0);
    double **Apk = allocate_2d(N,N,0);
    double **rkp1 = allocate_2d(N,N,0);
    
  
    for(int i = i1; i < i2; i++)
     {
      for(int j = i1; j < i2; j++)
    {
        uk[i][j] = 0;
        rk[i][j] = F[i][j];
        pk[i][j] = rk[i][j];    // # Start with uk = 0 --> r = b - Au = b            
    }   
     }   
   
   
 
    for(int k = 0; k < kmax; k++)
    {
        
        matvec(N, pk, Apk);
        
        double rTr = 0;
        double pTAp = 0;
        for(int i = i1; i < i2; i++)
        {
            for (int j = i1; j<i2;j++)
            {
            rTr += rk[i][j]*rk[i][j];
            pTAp += pk[i][j]*Apk[i][j]; 
            }
            
        }
         
        
        
        double alpha = rTr/pTAp;
        
        double rpTrp = 0;
        double max_res = 0;
        
        for(int i = i1; i < i2; i++)
        {
            for(int j = i1; j < i2; j++)
            {
                    uk[i][j] = uk[i][j] + alpha*pk[i][j];
                rkp1[i][j] = rk[i][j] - alpha*Apk[i][j];
                rpTrp += rkp1[i][j]*rkp1[i][j];
                max_res = fmax(fabs(rkp1[i][j]),max_res);
            }
        }
    
        
        itcount = k+1;
        
        //printf("%d error = %12.4e\n",k, max_res); 
        if (max_res < tol)
            break;
        
        
        double beta = rpTrp/rTr;
        
        for(int i = i1; i < i2; i++)
        {
         for(int j = i1; j < i2; j++)
            {
              pk[i][j] = rkp1[i][j] + beta*pk[i][j];   
            }
                  
        }
        
       
        for(int i = i1; i < i2; i++)
        {
         for (int j =i1; j<i2;j++)
               rk[i][j] = rkp1[i][j];  
        }
           
        rTr = rpTrp;
        
    }
    
    for(int i = i1; i < i2; i++)
    {
      for(int j = i1; j < i2; j++)
        {
          u[i][j] = uk[i][j];    
        }
             
    }
    
    
    
    free_2d(&uk,0);
    free_2d(&pk,1);
    free_2d(&rk,0);
    free_2d(&rkp1,0);
    free_2d(&Apk,0);

    
    
    return itcount;
    
}

int main(int argc, char** argv)
{

    int N;
    if (argc != 6)
    {
        N = 64;
        for(int i = 0; i < 4; i++)
            bc_type[i] = DIRICHLET;
    }        
    else
    {
        N = atoi(argv[1]);    
        for(int i = 0; i < 4; i++)
            bc_type[i] = atoi(argv[2+i]);
    }    
    
    // # Domain
    double a = 0, b = 2*M_PI;
    
    // # Numerical parameters
    double tol = 1e-13;
    int kmax = 10000;    
    int prt = 0;

    // # TODO : Implement a 2d solver using either cell-centered or node-centered. 
    
    double h = (b-a)/N;
    
    // # Arrays    
    double **u = allocate_2d(N,N,1);
    double **F = allocate_2d(N,N,0);
    double *x = allocate_1d(N,0);
    double *y = allocate_1d(N,0);
   
    // # Initialization
    for(int i = 0; i < N; i++)
    {
        x[i] = a + (i+0.5)*h;
        
        for(int j = 0; j < N; j++)
        {
          y[j] = a + (j+0.5)*h;  
            F[i][j] = pow(h,2)*upp_true(x[i],y[j]);
        }
    }
    
     
    
    //# Apply BCs
    
    
     if (bc_type[LEFT] == DIRICHLET)
    {
         for(int i = 0; i < N; i++)
        {
             F[0][i] -=2* u_true(a,y[i]);
        }
           
    }
    
    if (bc_type[LEFT] == NEUMANN)
    {
         for(int i = 0; i < N; i++)
        {
             F[0][i] -=h*un_true(LEFT,a, y[i]);
        }
           
    }
    
    if (bc_type[RIGHT] == DIRICHLET)
    {
        for (int i = 0; i< N; i++)
        {
            F[N-1][i] -=2*u_true(b,y[i]);
        }
    }
    
    if (bc_type[RIGHT] == NEUMANN)
    {
        for (int i = 0; i< N; i++)
        {
            F[N-1][i] -=h*un_true(RIGHT,b, y[i]);
        }
    }

    if (bc_type[TOP] == DIRICHLET)
    {
        for (int i = 0; i< N; i++)
        {
            F[i][N-1] -=2*u_true(b,x[i]);
        }
    }
        
    if (bc_type[TOP] == NEUMANN)
    {
        for (int i = 0; i< N; i++)
        {
            F[i][N-1] -=h*un_true(TOP,b, x[i]);
        }
    }
    
    if (bc_type[BOTTOM] == DIRICHLET)
    {
        for (int i = 0; i< N; i++)
        {
            F[i][0] -=2*u_true(a,x[i]);
        }
    }
        
    if (bc_type[BOTTOM] == NEUMANN)
    {
        for (int i = 0; i< N; i++)
        {
            F[i][0] -=h*un_true(BOTTOM,a, x[i]);
        }
    }
    
    // # Compute the solution using conjugate Gradient method
    int itcount = cg(N,F,u,tol,kmax,prt);
    
    double max_err = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j<N; j++)
        {
            max_err = fmax(fabs(u_true(x[i],y[j]) - u[i][j]),max_err);
            
        }
      
         
    }
 
     printf("max error = %12.4e\n", max_err); 
     
     FILE* fout = fopen("err_04.dat","w");   
     fwrite(&max_err,sizeof(double),1,fout);
     fwrite(&itcount,sizeof(int),1,fout);
     fclose(fout);
    
     
    FILE *file = fopen("cg_04.dat","wb");      
    fwrite(&N,sizeof(int),1,file);
    fwrite(bc_type,sizeof(int),4,file);
    fwrite(&a,sizeof(double),1,file);
    fwrite(&b,sizeof(double),1,file);
    fwrite(x,sizeof(double),N, file); 
    fwrite(y,sizeof(double),N, file); 
    
    for(int i = 0; i < N; i++)
        fwrite(u[i],sizeof(double),N, file); 
    fclose(file);
    
    
    free_2d(&u,1);
    free_2d(&F,0);
    free_1d(&x,0);
    free_1d(&y,0);

    
    return 0;
}    
