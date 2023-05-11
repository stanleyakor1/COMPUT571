
#include <stdio.h>

int main(int argc, char** argv)
{
    double x = 1.14; 
    
    printf("Value of x         (x) : %g\n",x);
    printf("Address of x      (&x) : %p\n",&x);
    printf("\n");
    
    double *y;
    y = &x;
    printf("Value of y         (y) : %p\n",y);
    printf("Value stored at y (*y) : %g\n",*y);
    
    return 0;
}
