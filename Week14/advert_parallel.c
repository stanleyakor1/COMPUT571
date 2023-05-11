#include "wpa1.c"

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

double * write(int N,int N_local,double **q_local,int rank)
{
    double *Q = allocate_1d(N*N,0);
    double *q1 = allocate_1d(N_local*N_local,0);
      for (int j = 0; j<N_local; j++)
    {
        for (int k = 0; k < N_local; k++)
        {
            q1[N_local*j+k] = q_local[j][k];
            
        }
        
    }
    
        

    MPI_Gather(q1, N_local*N_local, MPI_DOUBLE, Q, N_local*N_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    return Q;
}

int main(int argc, char** argv)
{
    FILE* fout = fopen("advect_parallel.dat","wb"); 
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int root = 0;
    int N;
    int ndim = 2;
    int dims[ndim];
    

    // # Number of processors in each direction.  Product should equal nprocs.
    dims[DIR_X] = atoi(argv[1]);
    dims[DIR_Y] = atoi(argv[2]);
        

    

    N = atoi(argv[3]);    
   

    if (rank == root)
        if (dims[DIR_X]*dims[DIR_Y] !=nprocs)
        {
            printf("dim[0]*dim[1] != nprocs\n");
            exit(0);
        }

    // # Domain
   
    double a = 0, b = 1;
    double Tfinal = 2;
    
    // # Numerical parameters
    int method_order = 2; 
    int prt = 0;
   
    MPI_Comm comm_cart;
    int periodicity[2] = {0,0};
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims,  periodicity, reorder, &comm_cart);

    // # Get coordinate of current rank
    int I, J, coords[ndim];
    {
        int maxdims = 2;
        MPI_Cart_coords(comm_cart,rank,maxdims,coords);
        I = coords[DIR_X];
        J = coords[DIR_Y];
        
    }
       
    // # Create domain  [0,P]x[0,1], where P=nprocs
    int P = nprocs;
   
    // # Set grid values
       
    // Assume Nx = Ny
    int Nx = N/dims[DIR_X];
    
    double dw = (b-a)/dims[DIR_X];
    double a_local_x = a + I*dw;
       
    double a_local_y = a + J*dw;

    double h = (b - a)/N;
            
    
   
     
    // # Set grid values and a true solution
    int N_local =Nx;
    // # Arrays    
    double **q_local = allocate_2d(N_local,N_local,1);
    double **uvel_local = allocate_2d(N_local+1,N_local,1);
    double **vvel_local = allocate_2d(N_local,N_local+1,1);
    
    double *xe_local= allocate_1d(N_local+1,1);
    double *xc_local= allocate_1d(N_local,1);
    
    // # Initialization
    
    int Imax = dims[DIR_X];
    int Jmax = dims[DIR_Y];
    
    for(int i = 0; i < N_local; i++)
    {

        xe_local[i] = a_local_x + i*h;
        xc_local[i] = a_local_x + (i+0.5)*h; 

    }
    
    xe_local[N_local] = a_local_x + N_local*h;
   
    
        // left ghost cells
        if (I == 0  )
        {
            
            xe_local[-1] = a_local_x - h;
            xc_local[-1] = a_local_x + (-1+0.5)*h;
                
        }
        else if(I == Imax-1)
        {
            xc_local[N_local] = a_local_x + (N_local+0.5)*h;
            xe_local[N_local+1] = b + h;
        }
        
     
       
     if (rank > 0)
    
    
        {
            int tag = 0;
            int sender = rank - 1;
            MPI_Recv(&xe_local[-1],1,MPI_DOUBLE,sender,tag,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
           
            
            tag = 1;
            int dest = rank - 1;
            MPI_Send(&xe_local[0],1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
          
           
            
            tag = 2;
            
            MPI_Recv(&xc_local[-1],1,MPI_DOUBLE,sender,tag,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         
            
            tag = 3;
           
            MPI_Send(&xc_local[0],1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
 
           
        }
        if (rank < nprocs - 1)

        {
            
            int tag = 0;
            int dest = rank +1;
            MPI_Send(&xe_local[N_local],1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
         
            
            tag = 1;
            int sender = rank + 1;
            MPI_Recv(&xe_local[N_local+1],1,MPI_DOUBLE,sender,tag,MPI_COMM_WORLD, MPI_STATUS_IGNORE);     
          
            
            tag = 2;
            MPI_Send(&xc_local[N_local-1],1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
          
            
            tag = 3;
         
            MPI_Recv(&xc_local[N_local],1,MPI_DOUBLE,sender,tag,MPI_COMM_WORLD, MPI_STATUS_IGNORE);     
         
            
        }
    

    
   
        // # Initialize q
        for(int i = -1; i < N_local+1; i++)
            for(int j = -1; j < N_local+1; j++)   // # initialize ghost cells
                q_local[i][j] = initial_condition(xc_local[i],xc_local[j]);

        // # Get u velocities at centes of edges of x faces
        double vvelx_local; // # not used at x-face
        for(int i = 0; i < N_local+1; i++)
            for(int j = 0; j < N_local; j++)            
                velocity(xe_local[i],xc_local[j],&uvel_local[i][j],&vvelx_local);

        // # Get v velocities at centesr of edges of y faces
        double uvely_local; // # not used at x-face
        for(int i = 0; i < N_local; i++)
            for(int j = 0; j < N_local+1; j++)            
                velocity(xc_local[i],xe_local[j],&uvely_local,&vvel_local[i][j]);

   
    // # Compute a time step
    double cfl = 0.45;
    double uvmax = 1;
    double dt_est = cfl*h/uvmax;
    int M = ceil(Tfinal/dt_est) + 1;
    double dt = Tfinal/M;
    
    double * Q = write(N, N_local,q_local,rank);
    double t = 0;
    
   
    double *qbig = NULL;
    if (rank == 0)
    {     
    fwrite(&N,sizeof(int),1,fout);    
    fwrite(&M,sizeof(int),1,fout);
    fwrite(&method_order,sizeof(int),1,fout);
    fwrite(&a,sizeof(double),1,fout);
    fwrite(&b,sizeof(double),1,fout);
    fwrite(&dt,sizeof(double),1,fout);
    
    // # write out initial condition
    
    fwrite(&t,sizeof(double),1,fout);
    qbig  = allocate_1d(N*N,0); 
    }
    
    // # Use parallel output routine above to gather solution to qbig on rank 0
    parallel_output(comm_cart,N,N_local,N_local,q_local,qbig);
    
    // # Write out the big array
    if (rank == 0)
    {
        fwrite(qbig,sizeof(double),N*N,fout);                         
        free_1d(&qbig,0);        
    }
     for(int n = 0; n < M; n++)
    {
        
        // # Update solution
        double cflmax = wpa_update(N_local, dt, h, uvel_local, vvel_local,  q_local, method_order);
        
        double MAXCFL;
        MPI_Reduce(&cflmax, &MAXCFL, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 
        
        double * Q = write(N, N_local,q_local,rank);
        
        if (rank == 0)
        {
          if (prt)
            printf("%12.4f %8.4f\n",t,MAXCFL);
        
            t += dt;
            fwrite(&t,sizeof(double),1,fout);
            qbig  = allocate_1d(N*N,0);
        }
        parallel_output(comm_cart,N,N_local,N_local,q_local,qbig);
        
        if (rank == 0)
        {
            fwrite(qbig,sizeof(double),N*N,fout);                         
            free_1d(&qbig,0);        
        }
    }     
    fclose(fout);
   
    //#endif
    free_2d(&q_local,1);
    free_2d(&uvel_local,1);
    free_2d(&vvel_local,1);
    free_1d(&xe_local,1);
    free_1d(&xc_local,1);
    
    MPI_Finalize();
    return 0;
}
