
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    int n = 5;
    
    size_t bytes = n*sizeof(double);
    
    double *x = (double*) malloc(bytes);  // cast when using with g++
    
    for(int i = 0; i < n; i++)
        x[i] = i*i;
    
    for(int i = 0; i < n; i++)
        printf("x[%d]  %6.0f\n",i,x[i]);
        
    free(x);
    return 0;
}
