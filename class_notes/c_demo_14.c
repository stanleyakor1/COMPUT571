
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    int n = 5;
    double *y;
    {
        double x[n];
        y = x;
    
        for(int i = 0; i < n; i++)
            x[i] = i;      
    }
    
    for(int i = 0; i < n; i++)
        printf("x[%d]  %6.1f\n",i,y[i]);
        
    return 0;
}
