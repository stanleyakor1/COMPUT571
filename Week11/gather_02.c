
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
    
    int root = 0;
    double * recv_buf;
    recv_buf = (double*) malloc(nprocs*sizeof(double));        
            
    int recv_count = 1;
    
    MPI_Gather(&r_pi,1,MPI_DOUBLE,recv_buf,recv_count,MPI_DOUBLE,root,MPI_COMM_WORLD);
    
    if (rank == root)
    {
        for(int p = 0; p < nprocs; p++)
            printf("Result from rank %d is %g\n",p,recv_buf[p]);        
    }        
    free(recv_buf);
    
    MPI_Finalize();
    return 0;
}
