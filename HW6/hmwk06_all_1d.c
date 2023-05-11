
enum 
{
    DIRICHLET = 0,
    NEUMANN,
};

enum
{
    LEFT=0,
    RIGHT
};

double* allocate_1d(int n, int m)
{
    double *mem = (double*) malloc((n + 2*m)*sizeof(double));
    return mem+m;
}

void free_1d(double **x, int m)
{
    free(*x-m);
    *x = NULL;
}

int parse_input(int argc, char** argv, int* N, int* bc_type)
{
    if (argc != 4)
    {
        *N = 64;
        bc_type[LEFT] = bc_type[RIGHT] = DIRICHLET;
    }
    else
    {
        *N = atoi(argv[1]);
        bc_type[LEFT] = atoi(argv[2]);
        bc_type[RIGHT] = atoi(argv[3]);
    }
    return 0;
}    


double u_true(double x)
{
    return (sqrt(M_PI)*x*erf(x) + exp(-x*x))/2.0;
}

double up_true(int iside, double x)
{
    double ux = sqrt(M_PI)*erf(x)/2.0;
    if (iside == LEFT)
        return -ux;
    else
        return ux;
}

double upp_true(double x)
{
    return exp(-x*x);
}
