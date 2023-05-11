
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "wpa.c"

double initial_condition(double x, double y)
{
    double y0 = 0.5;
    double r0 = 0.2;
    double r = fabs(y - y0);    
    if (r <= r0)
        return 0.25*(1 + cos(M_PI*r/r0));
    else
        return 0;
}

void velocity(double x, double y, double *u, double *v)
{
    *u = pow(sin(M_PI*x),2)*sin(2*M_PI*y);
    *v = -pow(sin(M_PI*y),2)*sin(2*M_PI*x);
}


int main(int argc, char** argv)
{
    int N = (argc == 2) ? atoi(argv[1]) : 64;
    
    // # Domain
    double a = 0, b = 1;
    double Tfinal = 2;
    
    int prt = 0;
    
    int method_order = 2;
    
    // # Arrays - Use cell centered mesh
    double **q = allocate_2d(N,N,1);
    double **uvel = allocate_2d(N+1,N,0);
    double **vvel = allocate_2d(N,N+1,0);

    // # Edges (xe) and centers (xc)
    double *xe = allocate_1d(N+1,1);
    double *xc = allocate_1d(N,1);
    
    // # Set up mesh;  include ghost cell values
    double h = (b-a)/N;    
    for(int i = -1; i < N+1; i++)
    {
        xe[i] = a + i*h;
        xc[i] = a + (i+0.5)*h;
    }
    xe[N+1] = b + h;

    // # Initialize q
    for(int i = -1; i < N+1; i++)
        for(int j = -1; j < N+1; j++)   // # initialize ghost cells
            q[i][j] = initial_condition(xc[i],xc[j]);
    
    // # Get u velocities at centes of edges of x faces
    double vvelx; // # not used at x-face
    for(int i = 0; i < N+1; i++)
        for(int j = 0; j < N; j++)            
            velocity(xe[i],xc[j],&uvel[i][j],&vvelx);

    // # Get v velocities at centesr of edges of y faces
    double uvely; // # not used at x-face
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N+1; j++)            
            velocity(xc[i],xe[j],&uvely,&vvel[i][j]);
    
    
    // # Compute a time step
    double cfl = 0.45;
    double uvmax = 1;
    double dt_est = cfl*h/uvmax;
    int M = ceil(Tfinal/dt_est) + 1;
    double dt = Tfinal/M;
    
    FILE* fout = fopen("advect.dat","wb");       
    fwrite(&N,sizeof(int),1,fout);    
    fwrite(&M,sizeof(int),1,fout);
    fwrite(&method_order,sizeof(int),1,fout);
    fwrite(&a,sizeof(double),1,fout);
    fwrite(&b,sizeof(double),1,fout);
    fwrite(&dt,sizeof(double),1,fout);
    
    // # write out initial condition
    double t = 0;
    fwrite(&t,sizeof(double),1,fout);
    for(int i = 0; i < N; i++)
        fwrite(q[i],sizeof(double),N, fout);             
    
    for(int n = 0; n < M; n++)
    {
        
        // # Update solution
        double cflmax = wpa_update(N, dt, h, uvel, vvel,  q, method_order);
        
        if (prt)
            printf("%12.4f %8.4f\n",t,cflmax);
        
        // # Write out current solution
        t += dt;
        fwrite(&t,sizeof(double),1,fout);
        for(int i = 0; i < N; i++)
            fwrite(q[i],sizeof(double),N, fout);             
    }        
    
    fclose(fout);
    
    free_2d(&q,1);
    free_2d(&uvel,0);
    free_2d(&vvel,0);
    free_1d(&xe,1);
    free_1d(&xc,1);
    
    return 0;
}    
