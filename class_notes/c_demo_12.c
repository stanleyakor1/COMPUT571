
#include <stdio.h>

int main(int argc, char** argv)
{
    double x = 1.2;
    {
        double x = 3.5;
    }
    printf("x = %f\n",x);
}
