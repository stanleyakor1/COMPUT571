
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <complex.h>
#include <cerf.h>

double* allocate_1d(int n, int m)
{
    /* Use malloc to allocate memory */
    double *y = (double*) malloc((n+2*m)*sizeof(double));
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
    /* TODO : Allocate memory for S and C */
    
     *C = allocate_1d(n+1,0);
     *S = allocate_1d(n+1,0);
     h = (bp-ap)/n;
    
    for(int i=0;i<=n;i++)
    {
       xi =(double) ap+i*h;
       *(*C+i)=(M_PI*(1-I))/8*(cerf(((1+I)/sqrt(2))*xi)+I*cerf(((1-I)/sqrt(2))*xi));
       *(*S+i)=(M_PI*(1+I))/8*(cerf(((1+I)/sqrt(2))*xi)-I*cerf(((1-I)/sqrt(2))*xi));
    }
    
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    int N = 1024;

    double a = -20.0, b = 20.0;
    
    int N_local = N/nprocs;

    /* TODO : Use "rank" to determine range of values [a_local,b_local] for this processor */
    
    double h = (b - a) / N;
    double a_local = a + rank * N_local * h;
    double b_local = a_local + N_local * h;

    
    double *C_local, *S_local;
    
    /* TODO : Compute range of values for this processor on [a_local, b_local] */
    
    eval_xy(a_local,b_local,N_local,&C_local,&S_local);
    
    if (rank == 0)
    {        
        /* Okay to define "automatic arrays" here, since they won't be referenced elsewhere */
        double C[N+1], S[N+1];
        
        /* Copy data from rank 0 into C and S*/
        for (int i = 0; i <= N_local; i++)
        {
            C[i] = C_local[i];
            S[i] = S_local[i];
        }
    
        
        /* TODO : Receive data from each of nprocs-1 processes into correct range in C, S */ 
        for (int i = 1; i < nprocs; i++) 
        {
         //TODO : Receive data from process i into correct range in C, S //
    
            MPI_Recv(&C[i * N_local], N_local+1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&S[i * N_local], N_local+1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }
        
        FILE *fout = fopen("prob8_output.dat","w");   
        
        /* TODO : Write out meta-data and arrays */
 
        fwrite(&N,sizeof(int),1,fout);
        fwrite(&a,sizeof(double),1,fout);
        fwrite(&b,sizeof(double),1,fout);
        fwrite(C,sizeof(double),N+1,fout);
        fwrite(S,sizeof(double),N+1,fout);
        fclose(fout);       
        
        
        /* ... */
         
    }
    else
    {

        /* TODO : Send data for C and S to rank 0 */
            MPI_Send(C_local, N_local+1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(S_local, N_local+1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }
    
    free_1d(&C_local);
    free_1d(&S_local);
    
    MPI_Finalize();
    return 0;
}
