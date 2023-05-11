
#include <stdio.h>

int main(int argc, char** argv)
{
    double x = 1.14; 
    double *y;
    y = &x;    /* Assign address to y */
    printf("Address of x : &x = %p\n",&x);
    printf("Value of y   :  y = %p\n",y);    
    return 0;
}
