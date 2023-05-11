
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    int n = 5;
    double *y;
    {
        double *x = (double*) malloc(n*sizeof(double));
        y = x;
    
        for(int i = 0; i < n; i++)
            x[i] = i;      
    }
    
    for(int i = 0; i < n; i++)
        printf("y[%d]  %6.1f\n",i,y[i]);
        
    free(y);
    return 0;
}
