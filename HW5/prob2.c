#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
double c2 = 0;
double c1 = -0.930763853398148;

double* allocate_1d(int n, int m)
{
    double *mem = (double*) malloc((n + 2*m)*sizeof(double));
    return &mem[m];
}

void free_1d(double **x, int m)
{
    free(&(*x)[-m]);
    *x = NULL;
}


double utrue(double x)
{
    double pi = M_PI;
    double utrue = (sqrt(pi)*x*erf(x) + exp(-x*x))/2.0 + c2*x + c1;
    return utrue;    
}

double rhs(double x)
{
    double upp = exp(-x*x);
    return upp;
}

void matvec(int N, double *u, double *L)
{
    
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

     if (rank > 0)
        {
            int tag = 0;
            int sender = rank - 1;
            MPI_Recv(&u[-1],1,MPI_DOUBLE,sender,tag,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            tag = 1;
            int dest = rank - 1;
            MPI_Send(&u[1],1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
            
           
        }
        if (rank < nprocs - 1)

        {
            
            int tag = 0;
            int dest = rank + 1;
            MPI_Send(&u[N-1],1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
            
            tag = 1;
            int sender = rank + 1;
            MPI_Recv(&u[N+1],1,MPI_DOUBLE,sender,tag,MPI_COMM_WORLD, MPI_STATUS_IGNORE);     
        }
    for(int i = 0; i <= N; i++)
    {
        L[i] = (u[i-1] - 2*u[i] + u[i+1]);
        
    }
}

void jacobi(int N, double *F, double *u, double tol, int kmax, int prt)
{
    
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    double *u_old = allocate_1d(N+1, 1);
    double *Lu = allocate_1d(N+1,0);
     double max_error = 0.0;
    
    double update;
    int i, j, iter;
    double error;
  
    for (i = 0; i <= N; i++) 
    {
        u_old[i] = u[i] = 0;
    }

    for (iter = 0; iter <= kmax; iter++) 
    { 

        if (rank == 0)
            {
                   u_old[-1] =  -u_old[1];
            }
            if (rank == nprocs -1)
            {
                   u_old[N+1] = - u_old[N-1];  
            }

        
        matvec(N,u_old,Lu);
        
        error = 0.0;
    
        for (i = 0; i < N+1; i++) 
            {
            update= -0.5*(F[i] - Lu[i]);
            u[i] = u_old[i]+update ;
            error = fmax(fabs(update), error);
           
            }
        MPI_Allreduce(&error, &max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (max_error < tol) 
        {
            break;
        }
        
          for (i = 0; i <N+1 ; i++) 
        {
            u_old[i]=u[i];
        }  
    }

    if (rank == 0 && prt == 1) 
    {
        printf("Number of iterations = %d\n", iter);
    }
    free_1d(&u_old,1);
    free_1d(&Lu,0);
}

int main(int argc, char **argv) 
{

    // Initialize MPI
     MPI_Init(&argc, &argv);
    
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Define problem parameters
     int N = atoi(argv[1]);
    double a = -1.0;
    double b = 1.0;
    double h = (b-a)/N;
    double tol = 1e-12;
    int kmax = 1e6;
    int prt = 1;

    // Define local variables for each process
   
    int N_local = N/nprocs;
    
    double a_local = a + rank * N_local * h;
    double b_local = a_local + N_local * h;
    
    double *F_local = allocate_1d(N_local+1,0);
    double *LU_local = allocate_1d(N_local+1,0);
    double *u_local = allocate_1d(N_local+1,1);
    double *utrue_local = allocate_1d(N_local+1,1);
    double h_local = (b_local - a_local)/N_local;
   
    for (int i =0; i<=N_local; i++)
    {
        double x = a_local + i*h_local;
        F_local[i] = pow(h_local,2)*rhs(x);
        utrue_local[i] = utrue(x); 
    } 
    jacobi(N_local, F_local, u_local, tol, kmax, prt);
    
    if (rank == 0)
    {
     u_local[-1] = -u_local[1];
    utrue_local[-1] = -utrue_local[1];
    }
        
    
    if (rank == nprocs -1)
    {
    u_local[N_local+1] = -u_local[N_local - 1];  
    utrue_local[N_local + 1] = -utrue_local[N_local - 1];
    }
        
    
    matvec(N_local,u_local,LU_local);
    
    double maxres = 0.0;
    double r = 0.0;
    for(int i = 0; i < N_local+1 ; i++)
    {
        r = fabs(F_local[i] - LU_local[i]);
        maxres = fmax(r,maxres);  
    }
    double maxerr = 0.0;
    double er = 0.0;
    for (int i =0 ; i<N_local+1; i++)
    {
        er = fabs(u_local[i] - utrue_local[i]);
        maxerr = fmax(er,maxerr);
    }
    
    MPI_Reduce(&er, &maxerr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&r, &maxres, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
     
      printf("Error norm = %12.4e\n", maxerr ); 
      printf("Maximum Residual = %12.4e\n", maxres );
      printf("-------------------------\n");
    }
        
     

    if (rank == 0)
    {        
        double U[N+1];
        for (int i = 0; i <= N_local; i++)
        {
            U[i] = u_local[i];
        }
        for (int i = 1; i < nprocs; i++) 
        {  
            MPI_Recv(&U[i * N_local], N_local+1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        FILE *fout = fopen("prob2.out","w");   
        fwrite(&N,sizeof(int),1,fout);
        fwrite(&a,sizeof(double),1,fout);
        fwrite(&b,sizeof(double),1,fout);
        fwrite(U,sizeof(double),N+1,fout);
        fclose(fout);       
         
    }
    else
    {
            MPI_Send(u_local, N_local+1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    
    free_1d(&F_local,0);
    free_1d(&u_local,1);
    free_1d(&LU_local,0);
    free_1d(&utrue_local,1);
    
   
    
    MPI_Finalize();

    return 0;
}
