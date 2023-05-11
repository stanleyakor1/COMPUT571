
#include <stdio.h>
#include <stdlib.h>

void f2(int n,double **y)
{
   *y = (double *)malloc(n*sizeof(double));
    printf("in f2 y=%p\n",*y);
    (*y)[0]=3.14159;
}

/* Function 2 goes here */

int main()
{
    double *y = NULL;
    int n = 10;
    f2(n,&y);
    printf("y[0] = %g\n",y[0]);
    printf("In main (after) y=%p\n",y);
    /* Call function 2 */
}
