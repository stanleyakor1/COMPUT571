
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* allocate_1d(int n)
{
    double *mem = (double*) malloc(n*sizeof(double));
    double *y = &mem[0];
    return y;
}

double  f(double x)
{
    return cos(M_PI*x);
}

int main(int argv, char** argc)
{
    int N = 128;
    
    double a = 0; 
    double b = 1; 
    
    double *x = allocate_1d(N+1);
    double *fx = allocate_1d(N+1);

    double h = (b-a)/N;
    
    for(int i = 0; i < N+1; i++)
    {
        x[i] = a + h*i;
        fx[i] = f(x[i]);
    }
    
    FILE *fout = fopen("array_01.dat","w");   
    fwrite(&N,sizeof(int),1, fout); 
    fwrite(&a,sizeof(double),1, fout); 
    fwrite(&b,sizeof(double),1, fout); 
    fwrite(&x[0],sizeof(double),N+1, fout); 
    fwrite(&fx[0],sizeof(double),N+1, fout);     
    fclose(fout);    
    
}
