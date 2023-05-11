
#include <stdio.h>     // Needed for file IO
#include <stdlib.h>    // Needed for 'atoi' and 'exit(0)'

int main(int argc, char** argv)
{
    if (argc < 4)
    {
        printf("Error : Include three magic ints at command line\n");
        exit(0);
    }
    int magic_int[3];
    
    for(int i = 0; i < 4; i++)
    {
        printf("%s\n",argv[i]);
    }

    for(int i = 0; i < 3; i++)
    {
        magic_int[i] = atoi(argv[i+1]);          
    }
    
    FILE *fout = fopen("magic_02.dat","w");   
    fwrite(&magic_int[0],sizeof(int),3, fout); 
    fclose(fout);

    return 0;
}
