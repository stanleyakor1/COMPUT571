
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
int main(int argc, char **argv) {

    // Initialize MPI
     MPI_Init(&argc, &argv);
    
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Define problem parameters
    int size  = atoi(argv[1]);
    
    int n = 100; // number of grid points
    double L = 1.0; // length of domain
    double dx = L / (n - 1); // grid spacing
    int n_local = n / size; // number of grid points per process
    double x_min = rank * n_local * dx; // minimum x value for this process
    double x_max = x_min + (n_local - 1) * dx; // maximum x value for this process
    
    
        double *u_local = (double*)malloc(n_local * sizeof(double));
    for (int i = 0; i < n_local; i++) {
        double x = x_min + i * dx;
        u_local[i] = sin(M_PI * x);
    }
    if (rank == 0) {
        u_local[0] = 0.0; // Dirichlet boundary condition at x = 0
    }
    if (rank == size - 1) {
        u_local[n_local - 1] = 0.0; // Dirichlet boundary condition at x = L
    }
    
        double tol = 1e-6; // convergence tolerance
    int max_iter = 100; // maximum number of iterations
    int iter = 0;
    double *u_old = (double*)malloc(n_local * sizeof(double));
    double *u_new = (double*)malloc(n_local * sizeof(double));
    while (iter < max_iter) {
        // save the old solution
        memcpy(u_old, u_local, n_local * sizeof(double));
        // exchange boundary data with neighboring processes
        if (rank > 0) {
            MPI_Sendrecv(u_local, 1, MPI_DOUBLE, rank - 1, 0, u_old - 1, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Sendrecv(u_local + n_local - 1, 1, MPI_DOUBLE, rank , 0, u_old + n_local, 1, MPI_DOUBLE, rank , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // compute the new solution
        for (int i = 0; i < n_local; i++) {
            double x = x_min + i * dx;
            u_new[i] = 0.5 * (u_old[i+1] + u_old[i-1] + dx*dx*sin(M_PI*x));
        }
        // check for convergence
        double local_error = 0.0;
        for (int i = 0; i < n_local; i++) {
            local_error += (u_new[i] - u_old[i]) * (u_new[i] - u_old[i]);
        }
        double global_error;
        MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        global_error = sqrt(global_error);
        if (global_error < tol) {
            break;
        }
        // update the solution
        memcpy(u_local, u_new, n_local * sizeof(double));
        iter++;
    }

        double *u = NULL;
    if (rank == 0) {
        u = (double*)malloc(n * sizeof(double));
    }
    MPI_Gather(u_local, n_local, MPI_DOUBLE, u, n_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        // print the solution
        for (int i = 0; i < n; i++) {
            double x = i * dx;
            printf("%f %f\n", x, u[i]);
        }
        free(u);
    }


    MPI_Finalize();

    
}
