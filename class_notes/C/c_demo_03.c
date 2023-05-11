
#include <stdio.h>

int main(int argc, char** argv)
{
    double x = 1.14;
    
    double y = x;    
    printf("y         = %g\n",y);
    
    x  = 5.4;
    printf("y (again) = %g\n",y);

    return 0;
}
