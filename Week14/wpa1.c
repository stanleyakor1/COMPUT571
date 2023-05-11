
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "advect_file1.c"

double wpa_update(int N,double dt, double h, double ** uvel, double **vvel,  
                  double** q, int order)
{
   
    double **apdq = allocate_2d(N+1,N,0);
    double **amdq = allocate_2d(N+1,N,0);

    double **bpdq = allocate_2d(N,N+1,0);
    double **bmdq = allocate_2d(N,N+1,0);
    
    
    // # X-sweep
    double uvmax = 0;
    for(int i = 0; i < N+1; i++)
        for(int j = 0; j < N; j++)
        {
            double wave = q[i][j] - q[i-1][j];
            double u = uvel[i][j];
            apdq[i][j] = fmax(u,0)*wave;
            amdq[i][j] = fmin(u,0)*wave;
            
            if (order == 2)
            {
                double cxx = 0.5*fabs(u)*(1 - dt*fabs(u)/h);
                apdq[i][j] -= cxx*wave;
                amdq[i][j] += cxx*wave;
            }
            uvmax = fmax(uvmax,fabs(u));
        }

    
    // # Y-sweep
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N+1; j++)
        {
            double wave = q[i][j] - q[i][j-1];
            double v = vvel[i][j];
            bpdq[i][j] = fmax(v,0)*wave;
            bmdq[i][j] = fmin(v,0)*wave;
            
            if (order == 2)
            {
                double cyy = 0.5*fabs(v)*(1 - dt*fabs(v)/h);
                bpdq[i][j] -= cyy*wave;
                bmdq[i][j] += cyy*wave;
            }
            uvmax = fmax(uvmax,fabs(v));
        }
    
    for(int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            q[i][j] -= dt/h*(apdq[i][j] + amdq[i+1][j] + bpdq[i][j] + bmdq[i][j+1]);        
    
    free_2d(&apdq,0);
    free_2d(&amdq,0);
    free_2d(&bpdq,0);
    free_2d(&bmdq,0);
    
    double cflmax = dt*uvmax/h;
    return cflmax;
}
