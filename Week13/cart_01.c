
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    int root = 0;
    
    // # Each processor stores a single integer
    int myval = rank + 10;
    
    int arr[nprocs];
    MPI_Gather(&myval,1,MPI_INTEGER,&arr[0],1,MPI_INTEGER,root,MPI_COMM_WORLD);
    if (rank == root)
    {
        printf("Before shift : \n");
        printf("%8s : ","Rank");
        for(int p = 0; p < nprocs; p++)
            printf("%5d",p);
        printf("\n");
        
        printf("%8s : ","Value");
        for(int p = 0; p < nprocs; p++)
            printf("%5d",arr[p]);
        printf("\n");
    }

    

    // #  Create a new communicator
    MPI_Comm comm_cart;  
    int ndim = 1;
    int dims[1] = {nprocs};
    int reorder = 0;
    
    
    int periodicity[1] = {0};
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims,  periodicity, reorder, &comm_cart);

    // # Shift information : 1: receive value from the left; -1
    int displ = 1;
    
    
    // #  source : This rank will receive data from SOURCE
    // #  dest   : This rank will send data to DESTINATION
    
    int dir = 0;
    int source, dest;
    MPI_Cart_shift(comm_cart, dir, displ, &source, &dest);
    
    int tag = 0;
    int myval_before = myval;
    // # Only a send buffer is required, no recv. buffer.
    MPI_Sendrecv_replace(&myval, 1, MPI_INTEGER, dest, tag, 
                                    source, tag, comm_cart, MPI_STATUS_IGNORE);

    if (rank == 0)
    {
        printf("\n");        
        
        if (displ== 1)
            printf("Values are shifted to the right (displ = 1). ",displ);
        else
            printf("Values are shifted to the left (displ = -1). ",displ);
        
        if (periodicity[0] == 1)
            printf("Domain is periodic\n\n");
        else
            printf("Domain is not periodic\n\n");
    }
    
    MPI_Gather(&myval,1,MPI_INTEGER,arr,1,MPI_INTEGER,root,MPI_COMM_WORLD);
    if (rank == root)
    {
        printf("After shift : \n");
        printf("%8s : ","Rank");
        for(int p = 0; p < nprocs; p++)
            printf("%5d",p);
        printf("\n");
        
        printf("%8s : ","Value");
        for(int p = 0; p < nprocs; p++)
            printf("%5d",arr[p]);
        printf("\n");
    }

    MPI_Finalize();
}
