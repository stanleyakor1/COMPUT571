
#include <stdio.h>

/* Function 1 goes here */
void f1(double *x)
{
    *x = 5.4;
}

int main()
{
    double x=1.3;
    printf("before : x = %g\n",x);
    f1(&x);
    printf("after : x = %g\n",x);
    
    return 0;
    /* Call function 1 */
}
