#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

double random_number()
{
  return (double) rand() / (double) RAND_MAX ;
}

void random_seed()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srand(clock() + rank);
}

int main(int argc, char** argv)
{
    int rank, nprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    random_seed();

    double x = random_number();
    printf("Rank %3d : Random number is %12.8f\n",rank,x);
    
    if (rank == 0)
    {         
          double xmax = x;
        
        for (int i = 1; i < nprocs; i++) 
        {
            double temp;
            MPI_Recv(&temp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            xmax = fmax(xmax, temp);
        }
        printf("Maximum value is %12.8f\n",xmax);       
    }
    else
    {
        MPI_Send(&x, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
}
