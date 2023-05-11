
#include <stdio.h>

int main(int argc, char** argv)
{
    double x = 1.2;
    printf("x = %g\n",x);
      
    double *y;
    y= &x;       
    *y = 3.5; 
    printf("x = %g\n",x);
    
    return 0;
}
