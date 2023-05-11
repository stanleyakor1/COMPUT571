
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <math.h>

void compute_stats(FILE *file, int n, double *mean, double *sum2)
{
    /* # TODO : Implement incremental algorithm for mean and sum2. 
       # Use only 'fread' to read 1 item at time from `file`.  
    */
    int i = 0;
    double dm, value;
    *mean = 0;
    *sum2 = 0;
    
    do {
        fread(&value, sizeof(double), 1, file);
        dm = value - *mean;
        *mean +=   dm/(double)(i + 1);
        *sum2 += dm * (value - *mean);
        i++;
    } while (i <n);
    
}


int main(int argc, char** argv)
{    
    FILE *file=fopen("X.dat","r");
    int N;
    fread(&N,sizeof(int),1,file);
    double mean, sum2;
    compute_stats(file,N,&mean, &sum2);
    fclose(file);

    printf("Mean  (C)       = %24.16f\n",mean);
    printf("Sum2  (C)       = %24.16f\n",sum2);
    printf("STD   (C)       = %24.16f\n",sqrt(sum2/N));
    
    return 0;
}
