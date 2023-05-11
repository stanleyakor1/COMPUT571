
#include <stdio.h>

int main(int argc, char** argv)
{
    double *y = NULL;    
    {
        double x = 11.5;
        y = &x;
    }        
    printf("y[0] = %f\n",y[0]);    
}
