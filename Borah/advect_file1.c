
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


double initial_condition(double x, double y)
{
    double y0 = 0.5;
    double r0 = 0.2;
    double r = fabs(y - y0);    
    if (r <= r0)
        return 0.25*(1 + cos(M_PI*r/r0));
    else
        return 0;
}

void velocity(double x, double y, double *u, double *v)
{
    *u = pow(sin(M_PI*x),2)*sin(2*M_PI*y);
    *v = -pow(sin(M_PI*y),2)*sin(2*M_PI*x);
}
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

double** allocate_2d(int n, int m, int mbc)
{
    int rows = n + 2*mbc;
    int cols = m + 2*mbc; 

    double   *qmem = malloc(rows*cols*sizeof(double));
    double **qrows = malloc(rows*sizeof(double*));

    for(int i = 0; i < rows; i++)
    {
        qrows[i] = &qmem[cols*i + mbc];
    }    
    return &qrows[mbc];
}

void free_2d(double ***q,int mbc)
{
    free(&(*q)[-mbc][-mbc]);
    free(&(*q)[-mbc]);
    *q = NULL;
}

enum
{
    DIR_X = 0,
    DIR_Y
};

enum
{
    LEFT=0,
    RIGHT,
    BOTTOM,
    TOP
};

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
