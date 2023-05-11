
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <math.h>

#include "cart_05_all.c"


void parallel_output(MPI_Comm comm_cart, int N, int Nx, int Ny, double **q, double *qbig)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int I,J;
    {
        int coords[2],ndims=2;
        MPI_Cart_coords(comm_cart,rank,ndims,coords);
        I = coords[DIR_X];
        J = coords[DIR_Y];
    }
    
    MPI_Comm sub_cart;
        
    // # Create a subgrid communicator.  Processors can only see those other procs that 
    // # remain in the subgrid.
    int remains[2] = {0,1}; // # Keep J in the communicator
    MPI_Cart_sub(comm_cart,remains,&sub_cart);
        
    // # Step 1 : Gather over J on each I subgrid.
    double *cols_local = (J == 0) ? allocate_1d(Nx*N,0) : NULL;
    for(int i = 0; i < Nx; i++)
        MPI_Gather(q[i],Ny,MPI_DOUBLE,&cols_local[N*i],Ny,MPI_DOUBLE,0,sub_cart); 
                
    // # Get a communicator for J=0
    MPI_Comm I_comm;
    MPI_Comm_split(MPI_COMM_WORLD,J,rank,&I_comm);
        
    // # Step 2 : Gather across I, at procs corresponding to J=0
    if (J == 0)
        MPI_Gather(cols_local,Nx*N,MPI_DOUBLE,qbig,Nx*N,MPI_DOUBLE,0,I_comm);      
    
    free_1d(&cols_local,0);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int root = 0;
    
    // #  Read input
    int ndims = 2;
    int dims[ndims];
    int N, Nx, Ny;
    if (argc == 4)
    {
        N = atoi(argv[1]);
        dims[DIR_X] = atoi(argv[2]);
        dims[DIR_Y] = atoi(argv[3]);
        
        // # Check that N is divisible by number of proces in each direction.
        Nx = N/dims[DIR_X];
        Ny = N/dims[DIR_Y];
        if (dims[DIR_X]*Nx != N || dims[DIR_Y]*Ny != N)
        {
            printf("N must be divisible by both dims[DIR_X] and dims[DIR_Y]\n");
            exit(0);        
        }        
    }
    else
        if (rank == 0)
        {
            printf("cart_05 <N> <dims 0> <dims 1>\n");
            exit(0);            
        }

    
    // # Create Cartesian communicator;  get coordinates of local processor.
    MPI_Comm comm_cart;
    int I, J;
    {
        int periodicity[2] = {0,0};
        int reorder = 0;
        MPI_Cart_create(MPI_COMM_WORLD, ndims, dims,  periodicity, reorder, &comm_cart);        
        
        int coords[2];
        MPI_Cart_coords(comm_cart,rank,ndims,coords);
        I = coords[DIR_X];
        J = coords[DIR_Y];
    }
    

    // # Set grid values and a true solution
    // # Note that data is contiguous over grid *columns*, i.e. constant i.  
    double a = 0, b = 1;
    double h = (b-a)/N;  
    double **q = allocate_2d(Nx,Ny,1);   // # Include ghost cells for illustration
    {
        for(int i = -1; i < Nx+1; i++)
            for(int j = -1; j < Ny+1; j++)
            {
                // # Row major ordering
                int i1 = Nx*I + i;
                int j1 = Ny*J + j;
                double x = a + i1*h;
                double y = a + j1*h;
                q[i][j] = qtrue(x,y); 
            }
    }
    
    // # Create big array for storing all of the data. 
    FILE *fout;
    double *qbig = NULL;
    if (rank == 0)
    {
        fout = fopen("cart_05.dat","wb");   
        fwrite(&N,sizeof(int),1,fout);
        fwrite(&a,sizeof(double),1,fout);
        fwrite(&b,sizeof(double),1,fout);        
        fwrite(&xshift,sizeof(double),1,fout);
        fwrite(&yshift,sizeof(double),1,fout);
        
        qbig  = allocate_1d(N*N,0);        
    }
    
    
    // # Use parallel output routine above to gather solution to qbig on rank 0
    parallel_output(comm_cart,N,Nx,Ny,q,qbig);
    
    // # Write out the big array
    if (rank == 0)
    {
        fwrite(qbig,sizeof(double),N*N,fout);             
        fclose(fout);            
        free_1d(&qbig,0);        
    }
        
    free_2d(&q,1);
    
    MPI_Finalize();
}
