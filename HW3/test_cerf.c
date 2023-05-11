
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <cerf.h>

int main(int argc, char** argv)
{
    complex double z0 = I*I;  /* should be -1 */
    printf("I*I = %20.16f + %20.16f\n",creal(z0),cimag(z0));
    
    /* Complex valued error function cerf */
    complex double z1 = cerf((1 + I)/sqrt(2.0));
    printf("z   = %20.16f + %20.16f\n",creal(z1),cimag(z1));

    return 0;
}    
