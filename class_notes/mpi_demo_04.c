
#include <mpi.h>
#include <stdio.h>    /* For IO */
#include <stdlib.h>   /* For malloc */

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Hardwire length of array */
    int N = 64;
    int Nlocal = N/nprocs;

    double *x;
    int tag = 0;
    if (rank == 0)
    {
        /* Allocate and initialize a big array for all data */
        x = malloc((N+1)*sizeof(double));
        for(int i = 0; i < N+1; i++)
        {
            x[i] = i;
        }
        for(int p = 1; p < nprocs; p++)
        {
            int dest = p;
            MPI_Send(&x[p*Nlocal],Nlocal+1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
        }
    }
    else
    {
        int source = 0;
        x = malloc((Nlocal+1)*sizeof(double));
        MPI_Recv(x,Nlocal+1,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    printf("Rank %d : (%2.0f,...,%2.0f)\n",rank,x[0],x[Nlocal]);
        
    free(x);
    
    MPI_Finalize();
}
