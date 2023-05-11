
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <math.h>

#include "cart_04_rhs.c"

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int root = 0;
    
    int N;
    if (argc == 2)
        N = atoi(argv[1]);
    else
        if (rank == 0)
        {
            printf("cart_04 <N>\n");
            exit(0);            
        }

    // # Domain
    int P = nprocs;
    double a = 0;
    double b = a + P;
    

    // # Create Cartesian communicator
    int ndim = 1;  // # 1 dimensional communicator (useful for 1d)
    int dims[ndim],periodicity[ndim];
    int reorder = 0;
    
    dims[DIR_X] = nprocs;
    periodicity[DIR_X] = 0;

    MPI_Comm comm_cart;
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims,  periodicity, reorder, &comm_cart);

    // # Get coordinate of current rank
    int I,coords[ndim];
    {
        int maxdims = 1;
        MPI_Cart_coords(comm_cart,rank,maxdims,coords);
        I = coords[DIR_X];
    }
    
    // # Create domain  [0,P]x[0,1], where P=nprocs
    double dw = (b-a)/nprocs;
    double a_local = a + I*dw;    
    double h = dw/N;  
    
    // # Set grid values and a true solution
    double *x = allocate_1d(N,1);
    double *y = allocate_1d(N,1);             
    double **q_soln = allocate_2d(N,N,1); 
    {
        for(int i = -1; i < N+1; i++)
        {
            x[i] = a_local + (i+0.5)*h;
            for(int j = -1; j < N+1; j++)
            {
                y[j] = a + (j + 0.5)*h;
                q_soln[i][j] = qtrue(x[i],y[j]);
            }
        }
    }
    
    double **F = allocate_2d(N,N,0);
    compute_RHS(N,q_soln,F,comm_cart);
    

    // # Compute a Laplacian at cell center.
    double maxerr = 0;
    double Ltrue = 1;
    double h2 = h*h;
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
        {
            maxerr = fmax(maxerr,fabs(F[i][j]/h2 -  Ltrue));        
        }
    
    if (rank == root)
    {
        printf("\n");
        printf("Maximum error :  %12.4e\n",maxerr);
        printf("\n");        
    }
    free_2d(&q_soln,1);
    free_2d(&F,0);
    free_1d(&x,1);
    free_1d(&y,1);
    
    MPI_Finalize();
}
