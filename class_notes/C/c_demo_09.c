
#include <stdio.h>

int main(int argc, char** argv)
{
    double x[3] = {1.,2.,3.}; 
    
    
    printf("x[0] = %f\n",x[0]);
    printf("x[1] = %f\n",x[1]);
    printf("x[2] = %f\n",x[2]);
    printf("\n");
    
    printf("Address of x %p\n",x);
    
    
    double *y = x; //x acts like a pointer
    printf("Value of y %p\n",y);
    y[1] = 15;
    
    printf("x[0] = %f\n",x[0]);
    printf("x[1] = %f\n",x[1]);
    printf("x[2] = %f\n",x[2]);
    
    return 0;
}
