
#include <stdio.h>

int main(int argc, char** argv)
{
    double *y;
    {
        double x;
        y = &x;
        
        x = 1.2;
    }    
    printf("*y = %24.16f\n",*y);        
}
