
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

enum
{
    LEFT=0,
    RIGHT
};

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    int root = 0;
    
    // #  Create a new communicator
    MPI_Comm comm_cart;  
    int ndim = 1;
    int dims = nprocs;
    int reorder = 0;        
    int periodicity[1] = {0};
    MPI_Cart_create(MPI_COMM_WORLD, ndim, &dims,  periodicity, reorder, &comm_cart);

    int maxdims = 2;
    int coords;
    MPI_Cart_coords(comm_cart, rank, maxdims, &coords);    
        
    // # Each processor stores a single integer
    int val = coords;  // # will be equal to rank in 1d case.
    
    int data[nprocs];
    MPI_Gather(&val,1,MPI_INTEGER,&data[0],1,MPI_INTEGER,root,MPI_COMM_WORLD);
    if (rank == root)
    {
        printf("Before shift : \n");
        printf("%8s : ","Rank");
        for(int p = 0; p < nprocs; p++)
            printf("%5d",p);
        printf("\n");
        
        printf("%8s : ","Value");
        for(int p = 0; p < nprocs; p++)
            printf("%5d",data[p]);
        printf("\n");
    }
    
    // # Shift information : 1: receive value from the left; -1 receive from right
    int displ = -1;
        
    int dir = 0;
    int nbr[2];
    MPI_Cart_shift(comm_cart, dir, displ, &nbr[LEFT], &nbr[RIGHT]);
    // # printf("rank : %d; source %d; dest %d\n",rank,nbr[LEFT],nbr[RIGHT]);
    
    int tag = 0;
    // int val_before = val;
    // # Only a send buffer is required, no recv. buffer.
    MPI_Sendrecv_replace(&val, 1, MPI_INTEGER, nbr[RIGHT], tag, 
                         nbr[LEFT], tag, comm_cart, MPI_STATUS_IGNORE);

    if (rank == 0)
    {
        printf("\n");        
        
        if (displ== 1)
            printf("Values are shifted to the right (displ = 1). %d\n",displ);
        else
            printf("Values are shifted to the left (displ = -1). %d\n",displ);
        
        if (periodicity[0] == 1)
            printf("Domain is periodic\n\n");
        else
            printf("Domain is not periodic\n\n");
    }
    
    MPI_Gather(&val,1,MPI_INTEGER,data,1,MPI_INTEGER,root,MPI_COMM_WORLD);
    if (rank == root)
    {
        printf("After shift : \n");
        printf("%8s : ","Rank");
        for(int p = 0; p < nprocs; p++)
            printf("%5d",p);
        printf("\n");
        
        printf("%8s : ","Value");
        for(int p = 0; p < nprocs; p++)
            printf("%5d",data[p]);
        printf("\n");
    }

    MPI_Finalize();
}
