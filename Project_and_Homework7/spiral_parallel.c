
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "cg_serial_fix.c"

double * write(int N,int N_local,double **q_local,int rank) //# gather file from all all processors
{
    double *Q = allocate_1d(N*N,0);
    double *q1 = allocate_1d(N_local*N_local,0);
      for (int j = 0; j<N_local; j++)
    {
        for (int k = 0; k < N_local; k++)
        {
            q1[N_local*j+k] = q_local[j][k];
            
        }
        
    }
    MPI_Gather(q1, N_local*N_local, MPI_DOUBLE, Q, N_local*N_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    free_1d(&q1,0);
    return Q;
}


int main(int argc, char** argv)
{
    FILE *fout = fopen("spiral_parallel.dat","w");
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int root = 0;
    int ndim = 2;
    int dims[ndim];
    if (argc > 3)
    {
        // # Number of processors in each direction.  Product should equal nprocs.
        dims[DIR_X] = atoi(argv[1]);
        dims[DIR_Y] = atoi(argv[2]);
        
    }
    else
    {
        printf("mpirun 4 spiral 2 2 <N> <nout> <Tfinal>\n");
        exit(0);  
    }
    
    int N, nout;
    double Tfinal;
    
    N = atoi(argv[3]);
    nout = atoi(argv[4]);  //# Number of time steps
    Tfinal = atof(argv[5]);
    
    
    for(int i = 0; i < 4; i++)
        bc_type[i] = NEUMANN;
    
    


    // Model parameters 
    double a_model = 0.75;
    double b_model = 0.01;
    double e_model = 0.02;

    // --------------------------- Numerical parameters -------------------------------

    int kmax = 1000;
    double tol = 1e-12;
    int prt = 0;
    // ---------------------------- Initialize solution -------------------------------
    
    if (rank == root)
        if (dims[DIR_X]*dims[DIR_Y] !=nprocs)
        {
            printf("dim[0]*dim[1] != nprocs\n");
            exit(0);
        }

    MPI_Comm comm_cart;
    int periodicity[2] = {0,0};
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims,  periodicity, reorder, &comm_cart);

    // # Get coordinate of current rank
    int I, J, coords[ndim];
    {
        int maxdims = 2;
        MPI_Cart_coords(comm_cart,rank,maxdims,coords);
        I = coords[DIR_X];
        J = coords[DIR_Y];
        
    }
   
    // # Set grid values
       
    // Assume Nx = Ny
    int Nx = N/dims[DIR_X];
    
    double L = 20;
    double a = -L, b = L;
    
    double dw = (b-a)/dims[DIR_X];
    double a_local_x = a + I*dw;
       
    double a_local_y = a + J*dw;

    double h = (b - a)/N;
    
    // # Set grid values and a true solution
    int N_local =Nx;
    // # Arrays 
    double **u = allocate_2d(N_local,N_local,1);
    double **v = allocate_2d(N_local,N_local,1);

    for(int i = -1; i < N_local+1; i++)
    {
        double x = a_local_x + (i+0.5)*h;
        for(int j = -1; j < N_local+1; j++)
        {
            double y = a_local_y + (j+0.5)*h;
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
    
    if (rank == 0)
        printf("dt = %f\n",dt);
    
    // # Time stepping
    double **up = allocate_2d(N_local,N_local,1);
    double **vp = allocate_2d(N_local,N_local,1);
    
    double t = 0;
    int Frame = 0;
    
    double *qbigU = NULL;
    double *qbigV = NULL;
    
    if (rank == 0)
    {
        printf("Frame %5d (step %5d)  t = %8.4f (itcount = %d)\n",Frame,0,t,0);       
        fwrite(&N,1,sizeof(int),fout);
        fwrite(&nout,1,sizeof(int),fout);
        fwrite(&a,1,sizeof(double),fout);
        fwrite(&b,1,sizeof(double),fout);
        fwrite(&t,1,sizeof(double),fout);    
        
        qbigU  = allocate_1d(N*N,0); 
        qbigV = allocate_1d(N*N,0); 
        
    }
    
    // # Use parallel output routine above to gather solution to qbig on rank 0
    parallel_output(comm_cart,N,N_local,N_local,u,qbigU);
    parallel_output(comm_cart,N,N_local,N_local,v,qbigV);
    
    // # Write out the big array
    if (rank == 0)
    {
        fwrite(qbigU,sizeof(double),N*N,fout); 
        fwrite(qbigV,sizeof(double),N*N,fout);       
        free_1d(&qbigU,0);
        free_1d(&qbigV,0);
    }

    
    double **F = allocate_2d(N_local,N_local,0);
    double lambda = h*h/dt;
    

    
    for(int n = 0; n < M+1; n++)
    {
        for(int i = 0; i < N_local; i++)
            for(int j = 0; j < N_local; j++)
            {
                double uij = u[i][j];
                double vij = v[i][j];
                double S = uij*(1-uij)*(uij - (vij+b_model)/a_model)/e_model;
                F[i][j] = -lambda*(u[i][j] + dt*S);
                up[i][j] = u[i][j];
            }
        
        int itcount = 0;
        itcount = cg(N_local,F,up,tol,kmax,prt,lambda, comm_cart);
        
       
        if (prt == 1)
            printf("Iteration count (CG) : %d\n",itcount);
        
        
        for(int i = 0; i < N_local; i++)
            for(int j = 0; j < N_local; j++)
                {
                    vp[i][j] = v[i][j] + dt*(up[i][j] - v[i][j]);
                }

        // # Write out current solution
        t += dt;
        if ((n+1)%M_inner == 0)
        {
            Frame++;
            
            if (rank == 0)
            {          
                printf("Frame %5d (step %5d)  t = %8.4f (itcount = %d)\n",Frame,n+1,t,itcount);  
                fwrite(&t,1,sizeof(double),fout); 
                qbigU  = allocate_1d(N*N,0); 
                qbigV = allocate_1d(N*N,0); 
            }  
            
            
            parallel_output(comm_cart,N,N_local,N_local,up,qbigU);
            parallel_output(comm_cart,N,N_local,N_local,vp,qbigV);

            // # Write out the big array
            if (rank == 0)
            {
                fwrite(qbigU,sizeof(double),N*N,fout); 
                fwrite(qbigV,sizeof(double),N*N,fout);       
                free_1d(&qbigU,0);
                free_1d(&qbigV,0);
            }
        
        }
       
        for(int i = 0; i < N_local; i++)
            for(int j = 0; j < N_local; j++)
            {
                u[i][j] = up[i][j];
                v[i][j] = vp[i][j];
            }
    }
    

    fclose(fout);
    free_2d(&u,1);
    free_2d(&v,1);
    free_2d(&up,1);
    free_2d(&vp,1);
    

    return 0;
}
