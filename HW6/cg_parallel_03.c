
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "hmwk06_all_1d.c"

int bc_type[2];

void matvec(int N, double *u, double *L)
{
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);    
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    int i1 = 0, i2 = N;   
    
    u[-1] = (rank == 0 && bc_type[LEFT] == NEUMANN)? u[0]:-u[0];
    u[N] = (rank == nprocs -1 && bc_type[RIGHT] == NEUMANN)? u[N-1]: -u[N-1]; 
    
    if (rank > 0)
    {
        int tag = 0;
        int sender = rank - 1;
        MPI_Sendrecv(&u[0], 1, MPI_DOUBLE, sender, tag, &u[-1], 1, MPI_DOUBLE, sender, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        

    }

    if (rank < nprocs - 1)
    {
        int tag = 0;
        int dest = rank + 1;
        MPI_Sendrecv(&u[N-1], 1, MPI_DOUBLE, dest, tag, &u[N], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }

    
    for(int i = i1; i < i2; i++) 
    {
        
            L[i] = u[i-1] - 2*u[i] + u[i+1];
    }

}

    
    
void cg(int N, double *F, double *u, double tol, int kmax, int prt,int *itcount)
{
    
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // # TODO : implement parallel CG
    int i1 =0, i2=N; 

    
    double *uk = allocate_1d(N,0);
    double *pk = allocate_1d(N,1);    
    double *rk = allocate_1d(N,0);
    double *Apk = allocate_1d(N,0);
    double *rkp1 = allocate_1d(N,0);
    
    
    for(int i = i1; i < i2; i++)
    {
        uk[i] = 0;
        rk[i] = F[i];
        pk[i] = rk[i];    // # Start with uk = 0 --> r = b - Au = b            
    }
    
    
    for(int i = i1; i < i2; i++)
    {
        uk[i] = 0;
        rk[i] = F[i];
        pk[i] = rk[i];    // # Start with uk = 0 --> r = b - Au = b            
    }
    
    for(int k = 0; k <kmax; k++)
    {        
        matvec(N,pk,Apk);
        
        double rTr = 0;
        double pTAp = 0;
        for(int i = i1; i < i2; i++)
        {
            rTr += rk[i]*rk[i];
            pTAp += pk[i]*Apk[i];
        }
        
        double rTr_global;
        double pTAp_global;
        
        
        MPI_Allreduce(&rTr, &rTr_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&pTAp, &pTAp_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
       
      
        double alpha = rTr_global/pTAp_global;
        
        double rpTrp = 0;
        double max_res = 0;
        for(int i = i1; i < i2; i++)
        {
            uk[i] = uk[i] + alpha*pk[i];
            rkp1[i] = rk[i] - alpha*Apk[i];
            rpTrp += rkp1[i]*rkp1[i];
            max_res = fmax(fabs(rkp1[i]),max_res);
        }
   
        double rpTrp_global;
        MPI_Allreduce(&rpTrp, &rpTrp_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double beta = rpTrp_global/rTr_global;
        
        *itcount = k+1;
        double MAXERROR;
        MPI_Allreduce(&max_res, &MAXERROR, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (prt)
            printf("%5d %12.4e\n",*itcount,MAXERROR);
        
        if (MAXERROR < tol)
        
            break;
                
        
        // # Update search directions
        for(int i = i1; i < i2; i++)
            pk[i] = rkp1[i] + beta*pk[i];
                
        // # Update residuals
        for(int i = i1; i < i2; i++)
            rk[i] = rkp1[i];  
    }
    for(int i = i1; i < i2; i++)
        u[i] = uk[i];     
    
    free_1d(&uk,0);
    free_1d(&pk,1);
    free_1d(&rk,0);
    free_1d(&Apk,0);
    free_1d(&rkp1,0);
        
}
    
    


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank,nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);    
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    int N;
    int count;
    int code = parse_input(argc,argv,&N,&bc_type[0]);
    if (code > 0)
        exit(1);
    
    // # Domain
    double a = 0.0, b = 1.0;
    
    // # Local values
    int N_local = N/nprocs;
    double h = (b - a)/N;
    double a_local = a + rank * N_local * h;
   
    double b_local = a_local + N_local * h;
    
   
    // # Numerical parameters
    int kmax;
    if (bc_type[LEFT] == NEUMANN && bc_type[RIGHT] == NEUMANN)
    {
        kmax = 30;
        printf("All Neumann boundary conditions may not converge. " \
               "Value of kmax is set to %d.\n\n",kmax);
    }
    else
        kmax = 10000;
        
    double tol = 1e-12;
    int prt = 0;
    
    // # TODO : Implemement parallel version of code in Problem 1 or Problem 2, above. 
    
    // # Set up right hand side F
    double h_local = (b_local - a_local)/N_local;
    
  
    
    // # Arrays    
    double *u_local = allocate_1d(N_local,1);
    double *F_local = allocate_1d(N_local,0);
    double *x_local = allocate_1d(N_local,0);
    
    // # Initialization
    for(int i = 0; i < N_local; i++)
    {
       x_local[i] = a_local + (i+0.5)*h_local;
       F_local[i] = h_local*h_local*upp_true(x_local[i]);       
    }

     
    
    if (rank == 0 && bc_type[LEFT] == DIRICHLET)
        F_local[0] -=2*u_true(a);
        
    
    if (rank == 0 && bc_type[LEFT] == NEUMANN)
    {
        F_local[0] -=h_local*up_true(LEFT,a);
    }
        
    
    if (rank == nprocs - 1 && bc_type[RIGHT] == DIRICHLET)
        F_local[N_local-1] -=2*u_true(b); 
        
     if (rank == nprocs -1 && bc_type[RIGHT] == NEUMANN)
        F_local[N_local-1] -=h_local*up_true(RIGHT,b);
 
     //# Compute the solution using CG
    cg(N_local,F_local,u_local,tol,kmax,prt, &count);
    
    double MAXERR = 0.0;
    double max_err = 0;
    for(int i = 0; i < N_local; i++)
    {
      x_local[i] = a_local + (i+0.5)*h_local;
      max_err = fmax(fabs(u_true(x_local[i]) - u_local[i]),max_err);
    }
        
    
    
    MPI_Reduce(&max_err, &MAXERR, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 
    
    if (rank == 0)
    {
    FILE* fout = fopen("err_03.dat","wb");   
    fwrite(&MAXERR,sizeof(double),1,fout);
    fwrite(&count,sizeof(int),1,fout);
    fclose(fout);   
    printf("itcount = %d\n",count);
    printf("\n");
    printf("N = %d\n",N);
     printf("Error (pde)     : %12.4e\n",MAXERR);
    }
              
    double *U = NULL, *x = NULL;
    if (rank == 0)
    {
        U = allocate_1d(N,0);
        x = allocate_1d(N,0);
    }
    MPI_Gather(&u_local[0], N_local, MPI_DOUBLE, U, N_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(x_local, N_local, MPI_DOUBLE, x, N_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        FILE *fout = fopen("cg_03.dat","w");
        fwrite(&N,sizeof(int),1,fout);
        fwrite(bc_type,sizeof(int),2,fout);
        fwrite(&a,sizeof(double),1,fout);
        fwrite(&b,sizeof(double),1,fout);
        fwrite(x,sizeof(double),N, fout);
        fwrite(U,sizeof(double),N, fout);
        fclose(fout);
    }
    free_1d(&F_local,0);
    free_1d(&u_local,1);
    free_1d(&x_local,0);
    free_1d(&U,0);
    free_1d(&x,0);
    
    MPI_Finalize();


    
    return 0;
}    
