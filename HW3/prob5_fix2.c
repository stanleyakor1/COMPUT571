
#include <stdio.h>  // second fix //

void change_value(double **y)
{
    
    **y = 11.5;
}

int main(int argc, char** argv)
{
    double x = 2.3;
    double *y = &x;
    
    printf("Old value : y[0] = %f\n",y[0]);
    change_value(&y);
    printf("New value : y[0] = %f\n",y[0]);    
}
