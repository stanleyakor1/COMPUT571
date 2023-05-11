#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define DIM 2

int main(int argc, char **argv) {
    int rank, size;
    int dims[DIM], periods[DIM], coords[DIM];
    int reorder = 0; // No reorder
    MPI_Comm comm_cart;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Define the dimensions of the 2D mesh
    dims[0] = atoi(argv[1]);; // number of rows
    dims[1] = atoi(argv[2]);; // number of columns
    
    if (rank == 0)
        {
            if (dims[0]*dims[1] != size)
            {
               printf("Usage : \n");
                printf("cart_03 <dim 0> <dim 1>\n");
                exit(0); 
            }
                        
        }
    // Set periodic boundary conditions (wraparound)
    periods[0] = 0; // row dimension is periodic
    periods[1] = 0; // column dimension is periodic

    // Create the Cartesian communicator
    MPI_Cart_create(MPI_COMM_WORLD, DIM, dims, periods, reorder, &comm_cart);

    // Get the coordinates of the current process in the 2D mesh
    MPI_Cart_coords(comm_cart, rank, DIM, coords);

    printf("Process %d: (%d,%d)\n", rank, coords[0], coords[1]);

    MPI_Finalize();
    return 0;
}
