
#include <mpi.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    int root = 0;
    
    // # Fill in q[-1] and q[1] with data from neighbors */
    MPI_Comm comm_cart;
    int ndim = 1;
    int dims[1] = {nprocs};
    int periodicity[1] = {0};
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims,  periodicity, reorder, &comm_cart);

    
    double qmem[3];    
    double *q = &qmem[1];

    double h = 1.0/nprocs;  
    
    int mycoords[1];
    MPI_Cart_get(comm_cart,ndim,dims,periodicity,mycoords);
    
    // # Fill in cell centered mesh locations using "my_coords"
    q[0] = (mycoords[0]+0.5)*h;
    
    double arr[nprocs];
    MPI_Gather(q,1,MPI_DOUBLE,arr,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
    if (rank == root)
    {
        printf("f(x) = x (cell centers) : \n");
        printf("%8s : ","Rank");
        for(int p = 0; p < nprocs; p++)
            printf("%10d",p);
        printf("\n");
        
        printf("%8s : ","f(x)");
        for(int p = 0; p < nprocs; p++)
            printf("%10.4f",arr[p]);
        printf("\n");
    }        

    // #  Fill ghost value q[-1]  (shift values right)
    int dir = 0;
    int disp = 1;
    int source, dest;
    int ierr = MPI_Cart_shift(comm_cart, dir, disp, &source, &dest);
    
    //# Need to check if we have valid source terms
    if (source == MPI_PROC_NULL)
    {
        //#we are at a left boundary
        q[-1] = -h/2;
    }

    int tag = 0;
    MPI_Sendrecv(&q[0], 1, MPI_DOUBLE, dest, tag, 
                 &q[-1], 1, MPI_DOUBLE, source, tag, 
                 comm_cart, MPI_STATUS_IGNORE);

    // #  Fill ghost value q[1]
    disp = -1;  
    MPI_Cart_shift(comm_cart, dir, disp, &source, &dest);
    
    if (source == MPI_PROC_NULL)
    {
        //#we are at a right boundary
        q[1] = 1+h/2;
    }


    tag = 1;
    ierr = MPI_Sendrecv(&q[0], 1, MPI_DOUBLE, dest, tag, 
                        &q[1], 1, MPI_DOUBLE, source, tag, 
                        comm_cart, MPI_STATUS_IGNORE);

    // # Compute a derivative using left and right values
    double deriv = (q[1] - q[-1])/(2*h);
    
    MPI_Gather(&deriv,1,MPI_DOUBLE,arr,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
    if (rank == root)
    {
        printf("\n");
        printf("Derivative f'(x) = 1 (cell-centers): \n");
        printf("%8s : ","Rank");
        for(int p = 0; p < nprocs; p++)
            printf("%10d",p);
        printf("\n");
        
        printf("%8s : ","f'(x)");
        for(int p = 0; p < nprocs; p++)
            printf("%10.4f",arr[p]);
        printf("\n");
    }
    

    MPI_Finalize();
}
