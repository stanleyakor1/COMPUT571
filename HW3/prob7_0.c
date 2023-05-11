#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* allocate_1d(int n, int m)
{
    /* Use malloc to allocate memory */
    double *mem = (double*) malloc((n+2*m)*sizeof(double));
    double *y = &mem[0];
    return y;
}


void free_1d(double **x, int m)
{
    /* Use free to free memory;  Set value of pointer to NULL after freeing memory.  */
    if (*x != NULL) 
    {
        free(*x);
        *x = NULL;
    }
    
}

void eval_xy(double ap, double bp, int n, double **C, double **S)
{    
    double xi, h;
    /* TODO : Allocate memory for S and C */
     *C = allocate_1d(n+1,0);
     *S = allocate_1d(n+1,0);
     h = (bp-ap)/n;
   
    /* TODO : Evaluate parametric equations C[i], S[i], i 0,1,...,n  over interval [ap,bp]*/
    for(int i=0;i<=n;i++)
    {
       xi =(double) ap+i*h;
       *(*C+i)=cos(xi);
       *(*S+i)=sin(xi);
    }
    
}

int main(void)
{
    int N = 10;
    double *Cptr;
    double *Sptr;
    
    eval_xy(0,2*M_PI,N,&Cptr,&Sptr);
    for(int i = 0; i<=N;i++){
        printf("%f : %f\n",Cptr[i],Sptr[i]);
    }
    
    free_1d(&Cptr,0);
    free_1d(&Sptr,0);
}
