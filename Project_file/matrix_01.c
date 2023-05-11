
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <math.h>

enum
{
    DIR_X = 0,
    DIR_Y
};

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int root = 0;
    
    int ndims = 2;
    int dims[ndims], N, M;
    if (argc == 3)
    {
        dims[DIR_X] = atoi(argv[1]);
        dims[DIR_Y] = atoi(argv[2]);
        if (rank == root && dims[DIR_X]*dims[DIR_Y] != nprocs)
        {
            printf("dims[DIR_X]*dims[DIR_Y] != nprocs\n");
            exit(0);
        }
        else
        {
            N = dims[0];
            M = dims[1];
            if (rank == 0)
                printf("Matrix is %d x %d\n",N,M);
        }
    }
    else
        if (rank == 0)
        {
            printf("matrix_01 <dims 0> <dims 1>\n");
            exit(0);            
        }
    

    // # Create Cartesian communicator
    int periodicity[2] = {0,0};
    int reorder = 0;
    
    MPI_Comm comm_cart;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims,  periodicity, reorder, &comm_cart);

    // # Get coordinate of current rank
    int coords[2];
    MPI_Cart_coords(comm_cart,rank,ndims,coords);
    int I = coords[DIR_X];
    int J = coords[DIR_Y];
    
    // # Each processor stores a single entry in the matrix
    double A = M*I + J;
    
    printf("Rank %2d : A[%d,%d] = %.0f\n",rank,I,J,A);
    
    // # Maximum column sum (1-norm)
    {
        MPI_Comm sub_cart;

        // # Reduce down each column;  store result in first row
        int remains[2] = {1,0};  // # Each i-index is sent to different subgrid
        MPI_Cart_sub(comm_cart,remains,&sub_cart);
        
        int RowID;
        MPI_Comm_rank(sub_cart,&RowID);
        
        int root_row = 0;  // # Result is stored in row 0
        double col_sum;
        double abs_Aij = fabs(A);
        MPI_Reduce(&abs_Aij, &col_sum, 1, MPI_DOUBLE, MPI_SUM, root_row, sub_cart); 
       
        // # Split communicators into rows
        MPI_Comm row_comm;
        MPI_Comm_split(MPI_COMM_WORLD, RowID, rank, &row_comm);
        
        if (RowID == 0)
        {
            int root_col = 0;
            double one_norm;
            MPI_Reduce(&col_sum, &one_norm, 1,MPI_DOUBLE,MPI_MAX,root_col,row_comm);        
            if (rank == root_col)
                printf("one-norm : %4.0f\n",one_norm);                    
        }
    }
        
    // # Maximum row sum (inf-norm)
    {
        MPI_Comm sub_cart;
        
        // # Reduce across each row;  store result in first column
        int remains[2] = {0,1};  // # j-indices all get sent to single subgrid communicator  
        MPI_Cart_sub(comm_cart,remains,&sub_cart);
        
        int ColID;
        MPI_Comm_rank(sub_cart,&ColID);
        
        int root_col = 0;
        double row_sum;
        double abs_Aij = fabs(A);
        MPI_Reduce(&abs_Aij,&row_sum, 1, MPI_DOUBLE,MPI_SUM,root_col,sub_cart); 
        
        // # Split communicator into columns   
        MPI_Comm column_comm;
        MPI_Comm_split(MPI_COMM_WORLD,ColID,rank,&column_comm);
        
        if (ColID == 0)
        {
            // # Only reduce down first column; store result in row 0.
            int root_row = 0;
            double inf_norm;
            MPI_Reduce(&row_sum, &inf_norm, 1,MPI_DOUBLE,MPI_MAX,root_row,column_comm);        
            if (rank == root_row)
                printf("inf-norm : %4.0f\n",inf_norm);                    
        }
    }

    MPI_Finalize();
}
