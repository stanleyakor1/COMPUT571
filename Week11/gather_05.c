
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <math.h>

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


void set_values(int n, double a, double b, double *q)
{
    double pi = M_PI;
    double dx = (b-a)/n;
    for(int i = 0; i < n+1; i++)
    {
        double x = a +  i*dx;
        q[i] = cos(x);        
    }    
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int N = 1 << 5;    // # 2^10

    double a = 0;
    double b = 2*M_PI;

    int N_local = N/nprocs;

    double dw = (b-a)/nprocs;
    double a_local = a + rank*dw;
    double b_local = a_local + dw;
    
    /* # Initialize data */
    double *q_local = allocate_1d(N_local+1,0);

    set_values(N_local,a_local, b_local, q_local);
    
    int root = 0;
    int *recv_counts, *displs;
    double *q;
    if (rank == root)
    {
        q = allocate_1d(N+1,0);
    
        recv_counts = (int*) malloc(nprocs*sizeof(int));
        displs = (int*) malloc(nprocs*sizeof(int));
    
        for(int p = 0; p < nprocs; p++)
        {
            recv_counts[p] = N_local+1;
            displs[p] = p*N_local;
        }
    }

    MPI_Gatherv(q_local,N_local+1,MPI_DOUBLE,q,recv_counts,displs,MPI_DOUBLE,root,MPI_COMM_WORLD);        
    
    if (rank == root)
    {
        FILE* fout = fopen("gather_05.out","w");
        fwrite(&N,1,sizeof(int),fout);
        fwrite(&a,1,sizeof(double),fout);
        fwrite(&b,1,sizeof(double),fout);

        fwrite(q,N+1,sizeof(double),fout);
        fclose(fout);        
        free_1d(&q,0);
    }

    MPI_Finalize();
    
    return 0;
}
