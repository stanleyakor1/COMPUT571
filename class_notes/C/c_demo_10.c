
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    int n = 5;
    size_t bytes = n*sizeof(int);
    double *x = (double*) malloc(bytes);
    
    for(int i = 0; i < n; i++)
    {
        x[i] = i;
    }
    
    printf("x[n-1] = %f\n",x[n-1]);
    
    free(x);
    return 0;
}
