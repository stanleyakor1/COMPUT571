
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "cg_serial_fix.c"

int main(int argc, char** argv)
{
    int N, nout;
    double Tfinal;
    if (argc == 4)
    {
        N = atoi(argv[1]);
        nout = atoi(argv[2]);     // # Number of time steps
        Tfinal = atof(argv[3]);
    }
    else
    {
        printf("spiral <N> <nout> <Tfinal>\n");
        exit(0);
    }

    for(int i = 0; i < 4; i++)
        bc_type[i] = NEUMANN;
    
    double L = 20;
    double a = -L, b = L;


    // Model parameters 
    double a_model = 0.75;
    double b_model = 0.01;
    double e_model = 0.02;

    // --------------------------- Numerical parameters -------------------------------
    double h = (b-a)/N;

    int kmax = 1000;
    double tol = 1e-12;
    int prt = 0;
    // ---------------------------- Initialize solution -------------------------------

    double **u = allocate_2d(N,N,1);
    double **v = allocate_2d(N,N,1);

    for(int i = -1; i < N+1; i++)
    {
        double x = a + (i+0.5)*h;
        for(int j = -1; j < N+1; j++)
        {
            double y = a + (j+0.5)*h;
            u[i][j] = (y > 0) ? 1 : 0;           
            v[i][j] = (x > 0) ? 1 : 0;
        }
    }

    // ----------------------------- Compute time step ---------------------------------
    // # Compute a stable time step
    // # 1.  Estimate a stable time step 'dt_stable'.   This may not evenly divide Tfinal. 
    // # 2.  Compute a minimum number M of time steps we need to take.
    // # 3.  Divide Tfinal by M to get get a dt that is guaranteed smaller than dt_est and  
    // #     satisfies M*dt = Tfinal.

        
    double dt_est = h/30;
    
    
    double dT = Tfinal/nout;
    int M_inner = ceil(dT/dt_est) + 1;   // # Compute M to guarantee we hit Tfinal
    double dt = dT/M_inner;
    int M = nout*M_inner;
    printf("dt = %f\n",dt);
    
    // # Time stepping
    double **up = allocate_2d(N,N,1);
    double **vp = allocate_2d(N,N,1);

    // #  Write out meta data 
    FILE *fout = fopen("spiral.dat","w");        
    fwrite(&N,1,sizeof(int),fout);
    fwrite(&nout,1,sizeof(int),fout);
    fwrite(&a,1,sizeof(double),fout);
    fwrite(&b,1,sizeof(double),fout);

    double t = 0;
    int Frame = 0;
    fwrite(&t,1,sizeof(double),fout);    
    printf("Frame %5d (step %5d)  t = %8.4f (itcount = %d)\n",Frame,0,t,0);    
    for(int i = 0; i < N; i++)
        fwrite(u[i],N,sizeof(double),fout);        
    
    for(int i = 0; i < N; i++)
        fwrite(v[i],N,sizeof(double),fout);        
    
    //M+1
    double **F = allocate_2d(N,N,0);
    double lambda = h*h/dt;
    for(int n = 0; n < M+1; n++)
    {
        for(int i = 0; i < N; i++)
            for(int j = 0; j < N; j++)
            {
                double uij = u[i][j];
                double vij = v[i][j];
                double S = uij*(1-uij)*(uij - (vij+b_model)/a_model)/e_model;
                F[i][j] = -lambda*(u[i][j] + dt*S);
                //printf("F[%d][%d] = %f\n",i,j,F[i][j]);
                up[i][j] = u[i][j];
               //printf("up[%d][%d] = %f, S = %f\n",i,j,up[i][j],S);
            }
        
        int itcount = 0;
        itcount = cg(N,F,up,lambda,tol,kmax,prt);
        if (prt == 1)
            printf("Iteration count (CG) : %d\n",itcount);
                
        for(int i = 0; i < N; i++)
            for(int j = 0; j < N; j++)
                {
                    vp[i][j] = v[i][j] + dt*(up[i][j] - v[i][j]);
                    //printf("up[%d][%d] = %f\n",i,j,up[i][j]);
                }

        // # Write out current solution
        t += dt;
        if ((n+1)%M_inner == 0)
        {
            Frame++;
            printf("Frame %5d (step %5d)  t = %8.4f (itcount = %d)\n",Frame,n+1,t,itcount);    
            
            fwrite(&t,sizeof(double),1,fout);
            for(int i = 0; i < N; i++)
                fwrite(up[i],sizeof(double),N, fout);             
        
            for(int i = 0; i < N; i++)
                fwrite(vp[i],sizeof(double),N, fout);  
        }
        
        for(int i = 0; i < N; i++)
            for(int j = 0; j < N; j++)
            {
                u[i][j] = up[i][j];
                v[i][j] = vp[i][j];
                //printf("up = %f, vp = %f, i = %d, j = %d\n",up[i][j],vp[i][j], i, j);
            }
    }
    fclose(fout);

    free_2d(&u,1);
    free_2d(&v,1);

    return 0;
}
