
#include <stdio.h>
#include <stdlib.h>

void f2(int n,double **y)
{
   double *mem = (double *)malloc(n*sizeof(double));
    *y = &mem[1];
    double *yptr = *y;
    yptr[-1] =15;
    yptr[0] =3.145;
}

/* Function 2 goes here */

int main()
{
    double *y = NULL;
    int n = 10;
    f2(n,&y);
    printf("y[0] = %g\n",y[0]);
    printf("y[-1] y=%g\n",y[-1]);
    /* Call function 2 */
}
