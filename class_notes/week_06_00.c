
#include <stdio.h>

double f0(double x)
{
    double value;
    value = x*x;
    
    return value;
}

int main()
{
    /* Call function 0 */
    double x = 1.2;

    double fx = f0(x);

    printf("f(x,y) = %g\n",fx);
    return 0;
}
