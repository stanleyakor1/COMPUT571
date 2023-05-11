
#include <mpi.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>

enum
{
    LEFT=0,
    RIGHT
};

enum
{
    DIR_X=0
};

double f(double x)
{
    return x*x/2;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    int root = 0;
    
    // # Create Cartesian communicator
    MPI_Comm comm_cart;
    int ndim = 1;  // # 1 dimensional communicator (useful for 1d)
    int dims[ndim],periodicity[ndim];
    int reorder = 0;
    
    dims[DIR_X] = nprocs;
    periodicity[DIR_X] = 0;
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims,  periodicity, reorder, &comm_cart);
    
    // # Get coordinate of current rank
    int I;
    {
        int maxdims = 2;
        MPI_Cart_coords(comm_cart, rank, maxdims, &I);            
    }
    
    // # construct local array  and index as : q[-1],q[0],q[1]    
    double qmem[3], *q;
    {
        q = &qmem[1];
    }

    // # Width of each processor "domain".
    double a,b,h,x;
    {
        a = 0;     
        b = 1;
        h = (b-a)/dims[DIR_X];
        
        // # Fill in cell centered mesh locations using coordinate I (not rank).
        x = a + (I + 0.5)*h;  // # Value at center of processor domain
        q[0] = f(x);
    }
    
    
    // # Write out data
    int grid_coords[nprocs];
    double x_grid[nprocs]; 
    
    MPI_Gather(&I,1,MPI_INTEGER,grid_coords,1,MPI_INTEGER,root,MPI_COMM_WORLD);
    MPI_Gather(&x,1,MPI_DOUBLE,x_grid,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
    if (rank == root)
    {
        printf("(I) (x) coordinates for each rank\n");
        for(int p = 0; p < nprocs; p++)        
            printf("%5d (%d) (%8.4f)\n",p,grid_coords[p],x_grid[p]);
        printf("\n");        
    }        

    int nbr[2];   // # Store left/right neighbors
        
    // #  Fill ghost value at the left : q[-1]
    {
        int dir = DIR_X;  // # only direction 0 in 1d
        int disp = 1;
        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &nbr[LEFT], &nbr[RIGHT]);

        // # Need to check if we have valid source rank from the left
        if (nbr[LEFT] == MPI_PROC_NULL)
            q[-1] = f(x-h);
    
        // # Send values to the right;  receive from the left
        int tag = 0;
        MPI_Sendrecv(&q[0], 1, MPI_DOUBLE, nbr[RIGHT], tag, 
                     &q[-1], 1, MPI_DOUBLE, nbr[LEFT], tag, 
                     comm_cart, MPI_STATUS_IGNORE);        
    }

    // #  Fill ghost value at the right : q[1]
    {
        int dir = DIR_X;
        int disp = -1;  
        MPI_Cart_shift(comm_cart, dir, disp, &nbr[RIGHT], &nbr[LEFT]);

        // # Need to check if we have valid source at the right
        if (nbr[RIGHT] == MPI_PROC_NULL)
            q[1] = f(x+h);
    
        // # Send values to the left;  receive from the right
        int tag = 1;
        int ierr = MPI_Sendrecv(&q[0], 1, MPI_DOUBLE, nbr[LEFT], tag, 
                                &q[1], 1, MPI_DOUBLE, nbr[RIGHT], tag, 
                                comm_cart, MPI_STATUS_IGNORE);        
    }

    // # Compute a derivative using left and right values
    double L = (q[1] - 2*q[0] + q[-1])/(h*h);
    
    // # Write out Laplacian
    double data[nprocs];
    MPI_Gather(&L,1,MPI_DOUBLE,data,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
    if (rank == root)
    {
        printf("\n");
        printf("Laplacian of f(x) = 0.5*x^2\n");
        for(int p = 0; p < nprocs; p++)
            printf("%5d %8.4f\n",p,data[p]);
        printf("\n");        
    }
    

    MPI_Finalize();
}
