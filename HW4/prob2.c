#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

void compute_stats(FILE *file, int n, double *mean, double *sum2)
{
    int i = 0;
    double dm, value;
    *mean = 0;
    *sum2 = 0;

    do {
        fread(&value, sizeof(double), 1, file);
        dm = value - *mean;
        *mean += dm/(i + 1);
        *sum2 += dm * (value - *mean);
        i++;
    } while (i < n);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    FILE *file=fopen("X.dat","r");
    int N;
    fread(&N,sizeof(int),1,file);

    int N_local = N/nprocs;
    double mean_local, sum2_local;

    fseek(file, rank * N_local * sizeof(double), SEEK_CUR);
    compute_stats(file, N_local, &mean_local, &sum2_local);
    fclose(file);

    MPI_Send(&mean_local, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&sum2_local, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

    double mean, sum2;
    mean = 0;
    sum2 = 0;

    if (rank == 0)
    {
        int NT = N_local;
        mean = mean_local;
        sum2 = sum2_local;
        for (int i = 1; i < nprocs; i++)
        {
            MPI_Recv(&mean_local, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&sum2_local, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            double delta = mean_local - mean;
            mean += delta * N_local / ((i + 1) *N_local);
            sum2 += sum2_local +  (pow(delta,2)* NT * N_local)/((i + 1) *N_local);
            NT +=N_local;
        }

        printf("Mean  (C)       = %24.16f\n", mean);
        printf("Sum2  (C)       = %24.16f\n", sum2);
        printf("STD   (C)       = %24.16f\n",sqrt(sum2/N));
    }

    MPI_Finalize();

    return 0;
}
