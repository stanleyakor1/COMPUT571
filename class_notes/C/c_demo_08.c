
#include <stdio.h>

int main(int argc, char** argv)
{
    double x[3]; 

    for(int i = 0; i < 3; i++)
    {
        x[i] = i*3.14159;
    }
    
    printf("x[0] = %g\n",x[0]);
    printf("x[1] = %g\n",x[1]);
    printf("x[2] = %g\n",x[2]);
    printf("\n");
    printf("x[3] (out of bounds!) %g\n",x[3]);
    
    return 0;
}
