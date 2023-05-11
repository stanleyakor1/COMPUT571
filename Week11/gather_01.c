
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <math.h>


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    double r_pi = rank*M_PI;    
    
    if (rank == 0)
    {
        int tag = 0;
        double *results = (double*) malloc(nprocs*sizeof(double));
        results[0] = r_pi;  // # assign rank 0 entry.
        for(int p = 1; p < nprocs; p++)
        {
            int source = p;            
            MPI_Recv(&results[p],1,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        for(int p = 0; p < nprocs; p++)
        {
            printf("Results on rank %d is %f\n",p,results[p]);
        }
        free(results);
    }
    else
    {
        int dest = 0;        
        int tag = 0;
        MPI_Send(&r_pi,1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
