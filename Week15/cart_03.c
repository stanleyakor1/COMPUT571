
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <math.h>

enum
{
    LEFT=0,
    RIGHT,
    DOWN,
    UP
};

enum
{
    DIR_X = 0,
    DIR_Y
};

double f(double x, double y)
{
    return (x*x + y*y)/4.0;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int root = 0;
    
    int ndim = 2;
    int dims[ndim];
    if (argc == 3)
    {
        // # Number of processors in each direction.  Product should equal nprocs.
        dims[DIR_X] = atoi(argv[1]);
        dims[DIR_Y] = atoi(argv[2]);
    }
    else
        if (rank == 0)
        {
            printf("Usage : \n");
            printf("cart_03 <dim 0> <dim 1>\n");
            exit(0);            
        }

    if (rank == root)
        if (nprocs != dims[DIR_X]*dims[DIR_Y])
        {
            printf("dim[0]*dim[1] != nprocs\n");
            exit(0);
        }
        
    // # Create Cartesian communicator
    MPI_Comm comm_cart;
    int periodicity[2] = {0,0};
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims,  periodicity, reorder, &comm_cart);
    
    // # Create a 3x3 block of memory
    double qmem[9], *qrows[3], **q;
    {
        for(int i = 0; i < 3; i++)
            qrows[i] = &qmem[3*i + 1];
        q = &qrows[1];        
    }
    
    // # Get coordinate of current rank
    int I, J,coords[2];
    {
        int maxdims = 2;
        MPI_Cart_coords(comm_cart,rank,maxdims,coords);

        I = coords[DIR_X];
        J = coords[DIR_Y];                
    }
    
    // # Set grid values
    double hx,hy,x,y;
    {
        double a = 0, b = 1;    // # Assume domain is square.
        hx = (b-a)/dims[DIR_X];   
        hy = (b-a)/dims[DIR_Y];
        x = a + (I+0.5)*hx;
        y = a + (J+0.5)*hy;
    
        // # Value
        q[0][0] = f(x,y);           
    }
    
    // # write out all 2d coordinates.
    int grid_coords[2*nprocs];
    double xy[2] = {x,y};
    double xy_grid[2*nprocs];
    
    MPI_Gather(coords,2,MPI_INTEGER,grid_coords,2,MPI_INTEGER,root,MPI_COMM_WORLD);
    MPI_Gather(xy,2,MPI_DOUBLE,xy_grid,2,MPI_DOUBLE,root,MPI_COMM_WORLD);
    if (rank == root)
    {
        printf("%8s     (%2s,%2s),   (%8s,%8s)\n","rank","I","J","x","y");
        for(int p = 0; p < nprocs; p++)
        {
            int I = grid_coords[2*p];
            int J = grid_coords[2*p+1];
            double x = xy_grid[2*p];
            double y = xy_grid[2*p + 1];
            printf("%8d     (%2d,%2d),   (%8.4f,%8.4f)\n",p,I,J,x,y);   
        }
        printf("\n");        
    }        

    // # Store source/dest neighbors
    int nbr[4];

    // # Get values at left : q[-1][0]
    {
        int dir = DIR_X;    
        int disp = 1;
        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &nbr[LEFT], &nbr[RIGHT]);

        if (nbr[LEFT] == MPI_PROC_NULL)
            q[-1][0] = f(x-hx,y);
    
        int tag = 0;
        MPI_Sendrecv(&q[0][0], 1, MPI_DOUBLE, nbr[RIGHT], tag, 
                     &q[-1][0], 1, MPI_DOUBLE, nbr[LEFT], tag, 
                     comm_cart, MPI_STATUS_IGNORE);  
        
    }

    // # Get values at the right : q[1][0]
    {
        int dir = DIR_X;
        int disp = -1;  
        MPI_Cart_shift(comm_cart, dir, disp, &nbr[RIGHT], &nbr[LEFT]);

        if (nbr[RIGHT] == MPI_PROC_NULL)
            q[1][0] = f(x+hx,y);
    
        int tag = 1;
        int ierr = MPI_Sendrecv(&q[0][0], 1, MPI_DOUBLE, nbr[LEFT], tag, 
                                &q[1][0], 1, MPI_DOUBLE, nbr[RIGHT], tag, 
                                comm_cart, MPI_STATUS_IGNORE);        
    }

    // #  Fill value at the bottom : q[0][-1]
    {
        int dir = DIR_Y;    
        int disp = 1;
        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &nbr[DOWN], &nbr[UP]);

        if (nbr[DOWN] == MPI_PROC_NULL)
            q[0][-1] = f(x,y-hy);

        int tag = 2;
        MPI_Sendrecv(&q[0][0], 1, MPI_DOUBLE, nbr[UP], tag, 
                     &q[0][-1], 1, MPI_DOUBLE, nbr[DOWN], tag, 
                     comm_cart, MPI_STATUS_IGNORE);        
    }

    // #  Fill value at the  top : q[0][1]
    {
        int dir = DIR_Y;
        int disp = -1;  
        MPI_Cart_shift(comm_cart, dir, disp, &nbr[UP], &nbr[DOWN]);
    
        if (nbr[UP] == MPI_PROC_NULL)
            q[0][1] = f(x,y+hy);

        int tag = 3;
        int ierr = MPI_Sendrecv(&q[0][0], 1, MPI_DOUBLE, nbr[DOWN], tag, 
                                &q[0][1], 1, MPI_DOUBLE, nbr[UP], tag, 
                                comm_cart, MPI_STATUS_IGNORE);        
    }

    // # Compute a Laplacian at cell center.
    double L = (q[-1][0] + q[1][0] - 2*q[0][0])/(hx*hx) + 
               (q[0][-1] + q[0][1] - 2*q[0][0])/(hy*hy);
    
    // # Write out the Laplacian
    double data[nprocs];
    MPI_Gather(&L,1,MPI_DOUBLE,data,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
    if (rank == root)
    {
        printf("\n");
        printf("Laplacian of f = 0.5*(x^2 + y*2) (cell-centers): \n");
        for(int p = 0; p < nprocs; p++)
            printf("%5d %8.4f\n",p,data[p]);
        printf("\n");        
    }
    MPI_Finalize();
}
