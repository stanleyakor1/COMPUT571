#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

double* allocate_1d(int nsize, int m)
{
    /* Use malloc to allocate memory */
    double *mem = (double*) malloc((nsize+2*m)*sizeof(double));
    double *y = &mem[m];
    return y;
}

void free_1d(double **x)
{
    /* Use free to free memory;  Set value of pointer to NULL after freeing memory.  */
    if (*x != NULL) 
    {
        free(*x);
        *x = NULL;
    }
}

void eval_xy(double ap, double bp, int n, double **C, double **S)
{    
    double xi, h;
    
    /* Allocate memory for S and C */
    *C = allocate_1d(n+1,0);
    *S = allocate_1d(n+1,0);
    
    h = (bp-ap)/n;
    for(int i=0;i<=n;i++)
    {
        xi = ap + i*h;
        (*C)[i] = cos(xi);
        (*S)[i] = sin(xi);
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    int N = 32;
    double a = 0, b = 2*M_PI;
    int N_local = N/nprocs;
    
    /* Use "rank" to determine range of values [a_local,b_local] for this processor */
    double a_local = a + rank*N_local*(b-a)/N;
    double b_local = a + (rank+1)*N_local*(b-a)/N;
    
    double *C_local, *S_local;
    
    /* Compute range of values for this processor on [a_local, b_local] */
    eval_xy(a_local, b_local, N_local, &C_local, &S_local);
    
    if (rank == 0)
    {
        /* Allocate arrays for C and S to store the results */
        double C[N+1], S[N+1];
        
        /* Receive data from each of nprocs-1 processes into correct range in C, S */
        for (int p = 1; p < nprocs; p++)
        {
            int offset = p*N_local;
            int count = N_local+1;
            MPI_Recv(&C[offset], count, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&S[offset], count, MPI_DOUBLE, p, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        /* Write out meta-data and arrays */
        FILE *fout = fopen("prob7_output.dat","w");   
        
        fwrite(&N,sizeof(int),1,fout);
        fwrite(&a,sizeof(double),1,fout);
        fwrite(&b,sizeof(double),1,fout);
        fwrite(C,sizeof(double),N+1,fout);
        fwrite(S,sizeof(double),N+1,fout);
        
        fclose(fout);          
    }
    else
    {

        /* TODO : Send data for C and S to rank 0 */
            MPI_Send(C_local, N_local+1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(S_local, N_local+1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
