
#include <stdio.h>
#include <mpi.h>

#include <string.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0)
    {
        char msg[100] = "Greetings from processor 0!";
        
        int dest = 1;
        int tag = 0;
        int len = strlen(msg)+1;
        printf("Sending to rank 1\n");
        MPI_Send(msg,len,MPI_CHAR,dest,tag,MPI_COMM_WORLD);
        printf("Done sending!\n");
    }
    else
    {
        MPI_Status status;
        int sender = 0;
        int tag = 0;
        char msg[100];
        printf("Processor %d is waiting to receive a message!\n",my_rank);
        MPI_Recv(msg,100,MPI_CHAR,sender,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        printf("Rank %d received message : \"%s\"\n",my_rank,msg);
    }
    
    MPI_Finalize();
}
