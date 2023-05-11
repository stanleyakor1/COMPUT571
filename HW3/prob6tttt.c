#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {
    int rank, size, max_val, my_val;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(rank); // Initialize random seed based on rank
    my_val = rand();

    printf("Process %d generated random value %d\n", rank, my_val);

    // Compute maximum value across all processes
    MPI_Reduce(&my_val, &max_val, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("The largest random number generated is %d\n", max_val);
    } else {
        MPI_Send(&my_val, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
