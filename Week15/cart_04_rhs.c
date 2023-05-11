
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <math.h>


#include "cart_04_all.c"

void comm(int N, double **q, MPI_Comm comm_cart)
{
    // # Store source/dest neighbors
    int nbr[2];
    double *sendbuf = allocate_1d(N,0);
    double *recvbuf = allocate_1d(N,0);
    
    
    // # Get values at left : q[-1][j]
    {
        int dir = DIR_X;    
        int disp = 1;
        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &nbr[LEFT], &nbr[RIGHT]);
        
        for(int j = 0; j < N; j++)
            sendbuf[j] = q[N-1][j];  // # Data sent to the right processor
        
        int sendtag = 0;
        int recvtag = 0;
        MPI_Sendrecv(sendbuf, N, MPI_DOUBLE, nbr[RIGHT], sendtag, 
                     recvbuf, N, MPI_DOUBLE, nbr[LEFT], recvtag, 
                     comm_cart, MPI_STATUS_IGNORE);         

        if (nbr[LEFT] != MPI_PROC_NULL)
            for(int j = 0; j < N; j++)
                q[-1][j] = recvbuf[j];  // # Data sent to the right processor                
    }

    // # Get values at the right : q[N][j]
    {
        int dir = DIR_X;
        int disp = -1;  
        MPI_Cart_shift(comm_cart, dir, disp, &nbr[RIGHT], &nbr[LEFT]);

        for(int j = 0; j < N; j++)
            sendbuf[j] = q[0][j];  // # Data sent to the left processor

        int sendtag = 1;
        int recvtag = 1;
        int ierr = MPI_Sendrecv(sendbuf, N, MPI_DOUBLE, nbr[LEFT],sendtag, 
                                recvbuf, N, MPI_DOUBLE, nbr[RIGHT], recvtag, 
                                comm_cart, MPI_STATUS_IGNORE);        
        
        if (nbr[RIGHT] != MPI_PROC_NULL)
            for(int j = 0; j < N; j++)
                q[N][j] = recvbuf[j];
    }       
}


void matvec(int N, double **q, double **L, MPI_Comm comm_cart)
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
    
    int nbr[2];
    
    // # Get values at left : q[-1][j]
    if (I == 0)  
        for(int j = 0; j < N; j++)
            q[-1][j] = -q[0][j];     // # Assume Dirichlet conditions

    
    // # Get values at right : q[N][j]
    if (I == Imax-1) 
        for(int j = 0; j < N; j++)
            q[N][j] = -q[N-1][j];     // # Assume Dirichlet conditions
        
    // # Fill value at the bottom and top : q[i][-1], q[i][N]
    // # Recall : We are assuming a 1 x P topology, so we don't have any neighbors 
    // # Above or below this rank.
    {
        for(int i = 0; i < N; i++)
        {
            q[i][-1] = -q[i][0];
            q[i][N] = -q[i][N-1];
        }
    }
        
    // # Communicate internal boundary data
    comm(N, q, comm_cart);
        
    // # Compute Laplacian (with homogenenous data)
    for(int i = 0; i < N; i++)    
        for(int j = 0; j < N; j++)
            L[i][j] = q[i-1][j] + q[i+1][j] + q[i][j-1] + q[i][j+1] - 4*q[i][j]; 
}
 

void compute_RHS(int N, double **q, double **F, MPI_Comm comm_cart)
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
    
    // # Compute right hand side as result of applying five-point Laplacian to 
    // # known solution.
    double **qh = allocate_2d(N,N,1);
    
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            qh[i][j] = q[i][j];
        
    matvec(N,qh,F,comm_cart);   
    
    // # Fix RHS with inhomogeneous data    
    for(int j = 0; j < N; j++)
    {
        if (I == 0)
            F[0][j] +=  (q[-1][j] + q[0][j]);
        if (I == Imax-1)
            F[N-1][j] += (q[N-1][j] + q[N][j]); 
    }

    for(int i = 0; i < N; i++)
    {
        F[i][0] += (q[i][0] + q[i][-1]);
        F[i][N-1] += (q[i][N-1] + q[i][N]);
    }
}

