
#include <stdio.h>
#include <stdlib.h>

double* allocate_1d(int n, int m)
{
    /* Use malloc to allocate memory */
    double *y = (double*) malloc((n+2*m)*sizeof(double));
    return y;
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

int main(int argv, char**argc)
{
    int N = 8;
    int m = 1;
    double a = 0;
    double b = 1;
    
    double *x = allocate_1d(N+1,m);
    
    double h = (b-a)/N;

    printf("%5s %8s\n","i","x[i]");
    printf("--------------\n");
    for(int i = -m; i < N+m+1; i++)
    {
        x[i] = a + i*h;
        printf("%5d %8.4f\n",i,x[i]);
       
    }
        
    free_1d(&x,0);

    
    return 0;
}
