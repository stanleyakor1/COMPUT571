
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "RD_all_2d.c"

int bc_type[4];

void matvec(int N, double **u, double **L, double lambda)
{
    for(int j = 0; j < N; j++)
    {
        u[-1][j] = (bc_type[LEFT] == DIRICHLET)  ? -u[0][j] : u[0][j];
        u[N][j]  = (bc_type[RIGHT] == DIRICHLET) ? -u[N-1][j] : u[N-1][j];
    }

    for(int i = 0; i < N; i++)
    {
        u[i][-1] = (bc_type[BOTTOM] == DIRICHLET) ? -u[i][0] : u[i][0];
        u[i][N] = (bc_type[TOP] == DIRICHLET) ? -u[i][N-1] : u[i][N-1];
    }
                    
    for(int i = 0; i < N; i++)    
        for(int j = 0; j < N; j++)
        {
            L[i][j] = u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] - 4*u[i][j]; 
            L[i][j] -= lambda*u[i][j];
        }
}

int cg(int N, double **F, double **u, double lambda, double tol, int kmax, int prt)
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
        
        matvec(N, pk, Apk,lambda);
        
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
