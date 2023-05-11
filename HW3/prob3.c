
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void free_1d(double **x, int m)
{
    /* Use free to free memory;  Set value of pointer to NULL after freeing memory.  */
    if (*x != NULL) 
    {
        free(*x+m);
        *x = NULL;
    }
    
}

double *linspace(int a,int b,int N)
{  
    double *y = (double*) malloc((N+1)*sizeof(double));
    double h = (double)(b-a)/N;
    
    for(int i=0;i<=N;i++)
    {
       y[i]=a+(double)i*h;
    }
    return y;
}


int main(int argv, char** argc)
{
    int N=100;
    
    /* Your code goes here*/ 
    double *y = linspace(-20,20,N);  
    double *x = linspace(-20,20,N);
    double r;
    double F[N+1][N+1];
    
    for(int i=0;i<=N;i++)
    {
        for(int j=0;j<=N;j++)
        {
            r=sqrt((x[i]-5)*(x[i]-5) + (y[j]+5)*(y[j]+5));
            if (r ==0)
                F[j][i]= 1;
            else
                F[j][i]= sin(r)/r;
        }
    }
    
    FILE *file = fopen("prob3_output.dat","w");
    fwrite(&N,sizeof(int),1, file);
    fwrite(&x[0],sizeof(double),N+1, file);
    fwrite(&y[0],sizeof(double),N+1, file);
    fwrite(&F, sizeof(double), (N+1)*(N+1), file);
    fclose(file);
    
    free_1d(&y,0);
    free_1d(&x,0);
    return 0;
}
