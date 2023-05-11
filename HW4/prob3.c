
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <mpi.h>

#include <math.h>
double* allocate_1d(int n, int m)
{
    /* Use malloc to allocate memory */
    double *y = (double*) malloc((n+2*m)*sizeof(double));
    return &y[m];
}


void free_1d(double **x, int m)
{
    /* Use free to free memory;  Set value of pointer to NULL after freeing memory.  */
    if (*x != NULL) 
    {
        free(*x+m);
        *x = NULL;
    }
    
}


double utrue(double x)
{
    return cos(x);
}

double upp(double x)
{
    return -cos(x);
}

/* Evaluate true solution */
void compute_solution(int n, double a, double b, double *u)
{
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    int i1 = 0;
    int i2 = (rank == nprocs-1) ? n+1 : n;
    
    double dx = (b-a)/n;
    for(int i = i1; i < i2; i++)
    {
        double x = a + dx*i;
        u[i] = utrue(x);
        
    }    
}


void compute_rhs(int n, double a, double b, double *f)
{
    double dx = (b-a)/n;
    for(int i = 0; i < n+1; i++)
    {
        double x = a + (dx*i);        
        f[i] = upp(x);        
    }
}

void apply_Laplacian(int n, double a, double b, double *u, double *L)
{
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* # TODO : Evaluate the Laplacian L.  You will need to exchange data at the boundary of
       # each processor region. */
    double h = (b - a)/n;
    int i2 = (rank == nprocs-1) ? n: n-1;
    int i1 = 0;
    
    if (rank > 0)
    {
        int tag = 0;
        int sender = rank - 1;
        MPI_Recv(&u[-1],1,MPI_DOUBLE,sender,tag,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        tag = 1;
        int dest = rank - 1;
        MPI_Send(&u[0],1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
    }
    if (rank < nprocs - 1)
    {
        int tag = 0;
        int dest = rank + 1;
        MPI_Send(&u[n-1],1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
      
        tag = 1;
        int sender = rank + 1;
        MPI_Recv(&u[n],1,MPI_DOUBLE,sender,tag,MPI_COMM_WORLD, MPI_STATUS_IGNORE);  
    }
    
    if (rank == 0)
    {
        u[-1] = 2*u[0] - u[1];
        
        
    }
    
    if (rank == nprocs -1)
    {
        u[n+1] = 2*u[n] - u[n - 1];
    }
    
   
    
    for (int j = i1; j <= i2; j++)
    {
        L[j] = (u[j - 1] - 2*u[j] + u[j+1])/pow(h,2);
        
    }
    
}


double compute_grid_error(int n, double a, double b, double *u, double *v)
{
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    // Compute local maximum error
    double local_max_error = 0.0;
    int i1 = (rank == 0) ? 1: 0;
    int i2 = (rank == nprocs - 1) ? n : n-1;
    
    for (int i = i1; i < i2; i++) {
        double error = fabs(u[i] - v[i]);
        if (error > local_max_error) {
            local_max_error = error;
        }
    }
    
    // Find global maximum error
    double global_max_error;
    MPI_Allreduce(&local_max_error, &global_max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    return global_max_error;
}


int main(int argc, char** argv)
{    
    MPI_Init(&argc, &argv);
    
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    if (argc < 2)
    {
        printf("User must supply N as a command line argument.\n");
        exit(1);
    }
    
    int N = atoi(argv[1]);
    double a = 0, b = 2*M_PI;
    int N_local = N/nprocs;
    
    double h = (b - a) / N;
    
    double a_local = a + rank * N_local * h;
    double b_local = a_local + N_local * h;
    
    FILE *file=fopen("prob3.dat","r");
    fread(&a,sizeof(double),1,file);
    fread(&b,sizeof(double),1,file);
    fclose(file);        
    
    /* # TODO : Determine interval [a_p, b_p] for this processor */
    double *u_local = allocate_1d(N_local+1,0);
    double *L_local = allocate_1d(N_local+1,0);
    double *v_local = allocate_1d(N_local+1,0);
    
    /* # TODO : Create array of u values (use 'utrue') */
    compute_solution(N_local, a_local, b_local, u_local);
    apply_Laplacian(N_local,a_local, b_local, u_local, L_local);
    compute_rhs(N_local, a_local,b_local,v_local);
    
    double g;
    g = compute_grid_error(N_local,a_local, b_local, L_local, v_local);
    double g_all; 
    MPI_Reduce(&g, &g_all, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0)
    {    
        /* # TODO  : Use MPI_Reduce to compute norm of truncation error for entire mesh */
        
        FILE *fout = fopen("prob3.out","w");
        fwrite(&g_all,sizeof(double),1,fout);
        fclose(fout);
    }
    free_1d(&u_local,0);
    free_1d(&L_local,0); 
    free_1d(&v_local,0);  
 
    MPI_Finalize();
}    
