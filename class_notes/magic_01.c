
#include <stdio.h>
#include <stdlib.h>    // Needed for 'atoi' and 'exit(0)'

int main(int argc, char** argv)
{
    if (argc <= 1)
    {
        printf("Error : Include magic int at command line\n");
        exit(0);
    }
    int magic_int = atoi(argv[1]);  

    FILE *fout = fopen("magic_01.dat","w");        
    fwrite(&magic_int,sizeof(int),1, fout); 
    fclose(fout);

    return 0;
}
