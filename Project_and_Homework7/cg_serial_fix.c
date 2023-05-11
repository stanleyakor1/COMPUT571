
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "RD_all_2d.c"

void parallel_output(MPI_Comm comm_cart, int N, int Nx, int Ny, double **q, double *qbig)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int I,J;
    {
        int coords[2],ndims=2;
        MPI_Cart_coords(comm_cart,rank,ndims,coords);
        I = coords[DIR_X];
        J = coords[DIR_Y];
    }
    
    MPI_Comm sub_cart;
        
    // # Create a subgrid communicator.  Processors can only see those other procs that 
    // # remain in the subgrid.
    int remains[2] = {0,1}; // # Keep J in the communicator
    MPI_Cart_sub(comm_cart,remains,&sub_cart);
        
    // # Step 1 : Gather over J on each I subgrid.
    double *cols_local = (J == 0) ? allocate_1d(Nx*N,0) : NULL;
    for(int i = 0; i < Nx; i++)
        MPI_Gather(q[i],Ny,MPI_DOUBLE,&cols_local[N*i],Ny,MPI_DOUBLE,0,sub_cart); 
                
    // # Get a communicator for J=0
    MPI_Comm I_comm;
    MPI_Comm_split(MPI_COMM_WORLD,J,rank,&I_comm);
        
    // # Step 2 : Gather across I, at procs corresponding to J=0
    if (J == 0)
        MPI_Gather(cols_local,Nx*N,MPI_DOUBLE,qbig,Nx*N,MPI_DOUBLE,0,I_comm);      
    
    free_1d(&cols_local,0);
}


void comm(int N, double **u, MPI_Comm comm_cart)
{
    int nbr[4];
    double *sendbuf = allocate_1d(N,0);
    double *recvbuf = allocate_1d(N,0);
    
    int maxdims = 2;
    int coords[maxdims];
    int periodicity[maxdims];
    int dims[maxdims];        
    
     int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    MPI_Cart_get(comm_cart,maxdims,dims,periodicity,coords);
    
     
    
   // Fill left ghost cell
    {      
       
        int tag = 0;
        int disp = 1;
        int dir = 0;
        for(int j = 0; j < N; j++)
            sendbuf[j] = u[N-1][j];

        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &nbr[LEFT], &nbr[RIGHT]);

        MPI_Sendrecv(sendbuf, N, MPI_DOUBLE, nbr[RIGHT], tag, 
                     recvbuf, N, MPI_DOUBLE, nbr[LEFT], tag, 
                     comm_cart, MPI_STATUS_IGNORE);

        if (nbr[LEFT] != MPI_PROC_NULL)
            for(int j = 0; j < N; j++)
                u[-1][j]  = recvbuf[j];       
    }     

   // Fill right ghost cell
    {    
        int tag = 1;
        
        int disp = -1;
        int dir = 0;
        for(int j = 0; j < N; j++)
            sendbuf[j] = u[0][j];

        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &nbr[RIGHT], &nbr[LEFT]);

        MPI_Sendrecv(sendbuf, N, MPI_DOUBLE, nbr[LEFT], tag, 
                     recvbuf, N, MPI_DOUBLE, nbr[RIGHT], tag, 
                     comm_cart, MPI_STATUS_IGNORE);


        if (nbr[RIGHT] != MPI_PROC_NULL)
            for(int j = 0; j < N; j++)
                u[N][j]  = recvbuf[j];

    }

  // fill bottom ghost cell
    {
        int tag = 2;
        
        int disp = 1;
        int dir = 1;
        for(int i = 0; i < N; i++)
            sendbuf[i] = u[i][N-1];

        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &nbr[BOTTOM], &nbr[TOP]);

        MPI_Sendrecv(sendbuf, N, MPI_DOUBLE, nbr[TOP], tag, 
                     recvbuf, N, MPI_DOUBLE, nbr[BOTTOM], tag, 
                     comm_cart, MPI_STATUS_IGNORE);

        if (nbr[BOTTOM] != MPI_PROC_NULL)
            for(int i = 0; i < N; i++)
                u[i][-1]  = recvbuf[i];
    }

    // Fill top ghost cell
    {  
        int tag = 3;
        int disp = -1;
        int dir = 1;
        for(int i = 0; i < N; i++)
            sendbuf[i] = u[i][0];

        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &nbr[TOP], &nbr[BOTTOM]);

        MPI_Sendrecv(sendbuf, N, MPI_DOUBLE, nbr[BOTTOM], tag, 
                     recvbuf, N, MPI_DOUBLE, nbr[TOP], tag, 
                     comm_cart, MPI_STATUS_IGNORE);

        if (nbr[TOP] != MPI_PROC_NULL)
            for(int i = 0; i < N; i++)
                u[i][N]  = recvbuf[i];                

    }
    free_1d(&sendbuf,0);
    free_1d(&recvbuf,0);
}



void matvec(int N, double **q, double **L,double lambda, MPI_Comm comm_cart)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // # Get current coordinate in comm_cart
    int maxdims = 1;
    int coords[maxdims];
    int periodicity[maxdims];
    int dims[maxdims];        
    
    MPI_Cart_get(comm_cart,maxdims,dims,periodicity,coords);
    int I = coords[DIR_X];    
    int Imax = dims[DIR_X];
    
    int J = coords[DIR_Y];    
    int Jmax = dims[DIR_Y];
    
    int nbr[4];
    
    // # Get values at left : q[-1][j]
     
   
    
    if (I == 0)  
        for(int j = 0; j < N; j++)
        {
          if (bc_type[LEFT] == DIRICHLET)
                q[-1][j] = -q[0][j];
            else
                q[-1][j] = q[0][j];
        }
            
    
    // # Get values at right : q[N][j]
    if (I == Imax-1) 
        for(int j = 0; j < N; j++)
        {
          if (bc_type[RIGHT] == DIRICHLET)
                q[N][j] = -q[N-1][j]; 
            else
                 q[N][j] = q[N-1][j];  
        }
    
      
            
    // # Fill value at the bottom and top : q[i][-1], q[i][N]
    if (J == 0)  
        for(int i = 0; i < N; i++)
        {
          if (bc_type[BOTTOM] == DIRICHLET )
                q[i][-1] = -q[i][0];     
            else
               q[i][-1] = q[i][0];  
        }
            
    if (J == Jmax-1)  
        for(int i = 0; i < N; i++)
        {
          if (bc_type[TOP] == DIRICHLET )
                q[i][N] = -q[i][N-1];     
            else
                q[i][N] = q[i][N-1];  
        }
            
        // # Communicate internal boundary data
       comm(N, q, comm_cart);          
   
        
    // # Compute Laplacian (with homogenenous data)
    for(int i = 0; i < N; i++)    
        for(int j = 0; j < N; j++)
            {
              L[i][j] = q[i-1][j] + q[i+1][j] + q[i][j-1] + q[i][j+1] - 4*(q[i][j]); 
              L[i][j] -= lambda*q[i][j];
            }
     
}
 

int cg(int N, double **F, double **u, double tol, int kmax, int prt,double lambda,MPI_Comm comm_cart)
{
    // # TODO : Implement at 2d CG solver
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // # Get current coordinate in comm_cart
  
    

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
        //printf("rk[%d][%d] = %f\n", i,j,F[i][j]);
        pk[i][j] = rk[i][j];    // # Start with uk = 0 --> r = b - Au = b            
    }   
     }   
   
   
 
    for(int k = 0; k < kmax; k++)
    {
        
        matvec(N, pk, Apk,lambda, comm_cart);
        
        double rTr = 0;
        double pTAp = 0;
        double alocal[2];
        
        for(int i = i1; i < i2; i++)
        {
            for (int j = i1; j<i2;j++)
            {
            rTr += rk[i][j]*rk[i][j];
            pTAp += pk[i][j]*Apk[i][j]; 
            }
            
        }
        alocal[0] = rTr;
        alocal[1] = pTAp;
        
        
        
        double a_1[2];
        MPI_Allreduce(alocal,a_1,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        
        
        double alpha = a_1[0]/a_1[1];
      
        
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
        
        double b_local;
        
        MPI_Allreduce(&rpTrp,&b_local,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    
        double error;
        MPI_Allreduce(&max_res,&error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 
        
        itcount = k+1;
        
        if (error < tol)
            break;
        
        
        double beta = b_local/a_1[0];
        
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
