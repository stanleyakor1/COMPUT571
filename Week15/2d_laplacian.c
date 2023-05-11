
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "hmwk07_all_2d.c"


void comm(int N, double **u, MPI_Comm comm_cart)
{
    int nbr[4];
    double *sendbuf = allocate_1d(N,0);
    double *recvbuf = allocate_1d(N,0);
    
    int tag = 0;

    int source, dest;
    
   // Fill left ghost cell
    {        
        int disp = 1;
        int dir = 0;
        for(int j = 0; j < N; j++)
            sendbuf[j] = u[N-1][j];

        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &nbr[LEFT], &nbr[RIGHT]);

        MPI_Sendrecv(sendbuf, N, MPI_DOUBLE, dest, tag, 
                     recvbuf, N, MPI_DOUBLE, source, tag, 
                     comm_cart, MPI_STATUS_IGNORE);

        if (nbr[LEFT] != MPI_PROC_NULL)
            for(int j = 0; j < N; j++)
                u[-1][j]  = recvbuf[j];       
    }     

   // Fill right ghost cell
    {        
        int disp = -1;
        int dir = 0;
        for(int j = 0; j < N; j++)
            sendbuf[j] = u[0][j];

        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &nbr[RIGHT], &nbr[LEFT]);

        MPI_Sendrecv(sendbuf, N, MPI_DOUBLE, dest, tag, 
                     recvbuf, N, MPI_DOUBLE, source, tag, 
                     comm_cart, MPI_STATUS_IGNORE);


        if (nbr[RIGHT] != MPI_PROC_NULL)
            for(int j = 0; j < N; j++)
                u[N][j]  = recvbuf[j];

    }

  // fill bottom ghost cell
    {
        int disp = 1;
        int dir = 1;
        for(int i = 0; i < N; i++)
            sendbuf[i] = u[i][N-1];

        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &nbr[BOTTOM], &nbr[TOP]);

        MPI_Sendrecv(sendbuf, N, MPI_DOUBLE, dest, tag, 
                     recvbuf, N, MPI_DOUBLE, source, tag, 
                     comm_cart, MPI_STATUS_IGNORE);

        if (nbr[BOTTOM] != MPI_PROC_NULL)
            for(int i = 0; i < N; i++)
                u[i][-1]  = recvbuf[i];
    }

    // Fill top ghost cell
    {        
        int disp = -1;
        int dir = 1;
        for(int i = 0; i < N; i++)
            sendbuf[i] = u[i][0];

        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &nbr[TOP], &nbr[BOTTOM]);

        MPI_Sendrecv(sendbuf, N, MPI_DOUBLE, dest, tag, 
                     recvbuf, N, MPI_DOUBLE, source, tag, 
                     comm_cart, MPI_STATUS_IGNORE);

        if (nbr[TOP] != MPI_PROC_NULL)
            for(int i = 0; i < N; i++)
                u[i][N]  = recvbuf[i];                

    }
    free_1d(&sendbuf,0);
    free_1d(&recvbuf,0);
}


void matvec(int N, double **q, double **L, MPI_Comm comm_cart)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // # Get current coordinate in comm_cart
    int maxdims = 1;
    int coords[maxdims];
    int periodicity[maxdims];
    int dims[maxdims];        
    
    MPI_Cart_get(comm_cart,maxdims,dims,periodicity,coords);
    int I = coords[DIR_X];    
    int Imax = dims[DIR_X];
    
    int J = coords[DIR_Y];    
    int Jmax = dims[DIR_Y];
    
    int nbr[4];
    
    // # Get values at left : q[-1][j]
    if (I == 0)  
        for(int j = 0; j < N; j++)
        {
          if (bc_type[LEFT] == DIRICHLET)
                q[-1][j] = -q[0][j];
            else
                q[-1][j] = q[0][j];
        }
            
    
    // # Get values at right : q[N][j]
    if (I == Imax-1) 
        for(int j = 0; j < N; j++)
        {
          if (bc_type[RIGHT] == DIRICHLET)
                q[N][j] = -q[N-1][j]; 
            else
                 q[N][j] = q[N-1][j];  
        }
            
        
        
    
    // # Fill value at the bottom and top : q[i][-1], q[i][N]
    if (J == 0)  
        for(int i = 0; i < N; i++)
        {
          if (bc_type[BOTTOM] == DIRICHLET)
                q[i][-1] = -q[i][0];     
            else
               q[i][-1] = q[i][0];  
        }
            
                
    
    if (J == Jmax-1)  
        for(int i = 0; i < N; i++)
        {
          if (bc_type[TOP] == DIRICHLET)
                q[i][N] = -q[i][N-1];     
            else
                q[i][N] = q[i][N-1];  
        }
            
                 
    
    // # Communicate internal boundary data
    comm(N, q, comm_cart);
        
    // # Compute Laplacian (with homogenenous data)
    for(int i = 0; i < N; i++)    
        for(int j = 0; j < N; j++)
            L[i][j] = q[i-1][j] + q[i+1][j] + q[i][j-1] + q[i][j+1] - 4*q[i][j]; 
}
 





int cg(int N, double **F, double **u, double tol, int kmax, int prt,MPI_Comm comm_cart)
{
    // # TODO : Implement at 2d CG solver
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // # Get current coordinate in comm_cart
    int maxdims = 1;
    int coords[maxdims];
    int periodicity[maxdims];
    int dims[maxdims];        
    
    MPI_Cart_get(comm_cart,maxdims,dims,periodicity,coords);
    int I = coords[DIR_X];    
    int Imax = dims[DIR_X];

    int itcount = 0;
    int i1 =0;
    int i2 = N;
    
    double **uk = allocate_2d(N,N,0);
    double **pk = allocate_2d(N,N,1);    
    double **rk = allocate_2d(N,N,0);
    double **Apk = allocate_2d(N,N,0);
    double **rkp1 = allocate_2d(N,N,0);
    
  
    for(int i = i1; i < i2; i++)
     {
      for(int j = i1; j < i2; j++)
    {
        uk[i][j] = 0;
        rk[i][j] = F[i][j];
        pk[i][j] = rk[i][j];    // # Start with uk = 0 --> r = b - Au = b            
    }   
     }   
   
   
 
    for(int k = 0; k < kmax; k++)
    {
        
        matvec(N, pk, Apk, comm_cart);
        
        double rTr = 0;
        double pTAp = 0;
        double alocal[2];
        
        for(int i = i1; i < i2; i++)
        {
            for (int j = i1; j<i2;j++)
            {
            rTr += rk[i][j]*rk[i][j];
            pTAp += pk[i][j]*Apk[i][j]; 
            }
            
        }
        alocal[0] = rTr;
        alocal[1] = pTAp;
        
        double a_1[2];
        MPI_Allreduce(alocal,a_1,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        
        
        double alpha = rTr/pTAp;
        
        double rpTrp = 0;
        double max_res = 0;
        
        for(int i = i1; i < i2; i++)
        {
            for(int j = i1; j < i2; j++)
            {
                uk[i][j] = uk[i][j] + alpha*pk[i][j];
                rkp1[i][j] = rk[i][j] - alpha*Apk[i][j];
                rpTrp += rkp1[i][j]*rkp1[i][j];
                max_res = fmax(fabs(rkp1[i][j]),max_res);
            }
        }
        
        double b_local;
        
        MPI_Allreduce(&rpTrp,&b_local,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    
        double error;
        MPI_Allreduce(&max_res,&error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 
        
        itcount = k+1;
        
        printf("%d error = %12.4e\n",k, max_res); 
        if (error < tol)
            break;
        
        
        double beta = b_local/a_1[0];
        
        for(int i = i1; i < i2; i++)
        {
         for(int j = i1; j < i2; j++)
            {
              pk[i][j] = rkp1[i][j] + beta*pk[i][j];   
            }
                  
        }
        
       
        for(int i = i1; i < i2; i++)
        {
         for (int j =i1; j<i2;j++)
               rk[i][j] = rkp1[i][j];  
        }
           
        rTr = rpTrp;
        
    }
    
    for(int i = i1; i < i2; i++)
    {
      for(int j = i1; j < i2; j++)
        {
          u[i][j] = uk[i][j];    
        }
             
    }
    
    
    
    free_2d(&uk,0);
    free_2d(&pk,1);
    free_2d(&rk,0);
    free_2d(&rkp1,0);
    free_2d(&Apk,0);

    
    
    return itcount;
    
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int root = 0;
    int N;
    int ndim = 2;
    int dims[ndim];
    if (argc > 3)
    {
        // # Number of processors in each direction.  Product should equal nprocs.
        dims[DIR_X] = atoi(argv[1]);
        dims[DIR_Y] = atoi(argv[2]);
        
    }
    
    if (argc != 7)
    {
        N = 64;
        for(int i = 0; i < 4; i++)
            bc_type[i] = DIRICHLET;
    }        
    else
    {
        N = atoi(argv[3]);    
        for(int i = 0; i < 4; i++)
            bc_type[i] = atoi(argv[4+i]);
    }    

    if (rank == root)
        if (dims[DIR_X]*dims[DIR_Y] !=nprocs)
        {
            printf("dim[0]*dim[1] != nprocs\n");
            exit(0);
        }

    // # Domain
    double a = 0, b = 2*M_PI;
    
    // # Numerical parameters
    double tol = 1e-13;
    int kmax = 10000;    
    int prt = 1;
   
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
        
        //printf("Process %d: (%d,%d)\n", rank, coords[0], coords[1]);
    }
    
    // # Create domain  [0,P]x[0,1], where P=nprocs
    int P = nprocs;
   
    // # Set grid values
   
        double dw = (b-a)/P;
        double a_local_x = a + J*dw*dims[DIR_X];
        double b_local_x = a_local_x + dw*dims[DIR_X];
        double a_local_y = a + I*dw*N;
        double b_local_y = a_local_y + dw*N;

        //printf("a_local = %f, b_local = %f\n", a_local_x, b_local_x);

        double h_x = (b_local_x - a_local_x)/dims[DIR_X];
        double h_y = (b_local_y - a_local_y)/dims[DIR_Y];    
    
    
     
    // # Set grid values and a true solution
    int N_local =N/dims[DIR_X];
    // # Arrays    
    double **u = allocate_2d(N_local,N_local,1);
    double **F = allocate_2d(N_local,N_local,0);
    double *x = allocate_1d(N_local,0);
    double *y = allocate_1d(N_local,0);
    
    // # Initialization
    for(int i = 0; i < N_local; i++)
    {
        x[i] = a_local_x + (i+0.5)*h_x;
        for(int j = 0; j < N_local; j++)
        {  
            y[i] = a_local_y + (i+0.5)*h_y;
            F[i][j] = pow(h_x,2)*upp_true(x[i],y[j]);
        }
    }
    
    int Imax = dims[DIR_X];
    int Jmax = dims[DIR_Y];
    
    for(int i = 0; i < N_local; i++)
    {
        if (I == 0)
            {
            if (bc_type[LEFT] == DIRICHLET)
                F[0][i] -=  2* u_true(a_local_x,y[i]);
            else
                F[0][i] -=h_x*un_true(LEFT,a_local_x, y[i]);
            }
        
        if (I == Imax-1)
        {
           if (bc_type[RIGHT] == DIRICHLET)
                F[N_local-1][i] -=2*u_true(b_local_x,y[i]);
            else
                F[N_local-1][i] -=h_x*un_true(RIGHT,b_local_x, y[i]); 
        }
            
        
        if (J == 0)
        {
            if (bc_type[BOTTOM] == DIRICHLET)
                F[i][0] -=2*u_true(a_local_x,x[i]);
            else
                F[i][0] -=h_x*un_true(BOTTOM,a_local_x, x[i]);
            
        }
            
        
        if (J == Jmax-1)
        {
          if (bc_type[TOP] == DIRICHLET)
                F[i][N_local-1] -=2*u_true(b_local_x,x[i]);
            else
                F[i][N_local-1] -=h_x*un_true(TOP,b_local_x, x[i]);   
        }
              
    }
    
    int itcount = cg(N_local,F,u,tol,kmax,prt,comm_cart);
    //printf("itcount = %d\n", itcount);
    
    free_1d(&x,0);
    free_1d(&y,0);
    free_2d(&F,0);
    free_2d(&u,1);
    MPI_Finalize();
    return 0;
}
