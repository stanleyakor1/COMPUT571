
#include <mpi.h>
#include <stdio.h>    /* For IO */

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Hardwire length of array */
    int tag = 0;
    int data;
    if (rank == 0)
    {
        for(int p = 1; p < nprocs; p++)
        {
            int dest = p;
            int data = p;
            MPI_Send(&data,1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
        }
    }
    else
    {
        int source = 0;        
        MPI_Recv(&data,1,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    printf("Rank %d : (%d)\n",rank,data);
        
    MPI_Finalize();
}
