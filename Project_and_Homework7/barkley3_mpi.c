/*
 * The the 2D Barkley model
 * Written by:	Michael Chiwere
 * 		Boise State University
 * 		03/11/2021
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define ID_2D(i,j,nx) ((j)*(nx+2)+(i))

//** Function interfaces **

// getArgs_mpi reads input parameters Nfrom command line

void getArgs(int *N, int *L, double *a, double *b, double *eps, int *tfinal, int argc, char *argv[], int irank, MPI_Comm comm);

double f(double u, double v, double a , double b, double eps);

double g(double u, double v);

void write_results(double *u, double *v, int N_loc, int nproc, int irank, int N, int L, double h);

//** Main function **

int main(int argc, char *argv[])
{
	int N, L, i, j, id, id_ghost,tfinal;
	int id_lft, id_rgt, id_bot, id_top;
	double  eps, a, b;
        double t=0.0;

	int irank, nproc;

	//declare variables for tracking requests
 	MPI_Request request_send_rgt;
  	MPI_Request request_send_lft;
        MPI_Request request_send_bot;
        MPI_Request request_send_top;

      MPI_Request request_recv_rgt;
      MPI_Request request_recv_lft;
      MPI_Request request_recv_bot;
     MPI_Request request_recv_top; 
        
	//initialize MPI
  	MPI_Init(&argc,&argv); 
  	MPI_Comm_rank(MPI_COMM_WORLD,&irank);
  	MPI_Comm_size(MPI_COMM_WORLD,&nproc);

  //compute row and column index in the rank grid
  	int q = (int)sqrt(nproc); //assume nproc is a square number
  	int rank_row = irank/q;
  	int rank_col = irank%q;

	//read command line arguments
	getArgs( &N, &L, &a, &b,  &eps, &tfinal,  argc, argv, irank, MPI_COMM_WORLD);

	// compute the grid spacing
	
	double h =(2.0 *L)/(N-1);
	double k = 0.0005;
	// compute local problem size
	
	int N_loc = N/q;

	// allocate arrays
	
	double *u       = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  	double *unew   = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  	double *vnew = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  	double *v       = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  	double *x       = malloc((N_loc+2)*sizeof(double));
  	double *y       = malloc((N_loc+2)*sizeof(double));

	  //allocate buffers
	
       double *send_buffer_rgt= malloc(N_loc*sizeof(double));
	 double *send_buffer_lft = malloc(N_loc*sizeof(double));
 	 double *send_buffer_bot = malloc(N_loc*sizeof(double));
	 double *send_buffer_top = malloc(N_loc*sizeof(double));

 	double *recv_buffer_rgt = malloc(N_loc*sizeof(double));
 	double *recv_buffer_lft = malloc(N_loc*sizeof(double));
 	double *recv_buffer_bot = malloc(N_loc*sizeof(double));
       double *recv_buffer_top = malloc(N_loc*sizeof(double));

	

	// initialize array u
	
	 for(i=1;i<N_loc+1; i++){
                x[i] = -L + (i-1)*h + rank_col*N_loc*h;
                y[i] = -L + (i-1)*h + rank_row*N_loc*h;
        }
	
	 // initializing u and v  arrays

	  for(j=1; j<N_loc+1; j++){
                for(i=1; i<N_loc+1; i++){
                        id = ID_2D(i,j,N_loc);
                        if(y[j] >= 0)
                        {
                                u[id] = 1;
                        }
                        else
                        {
                                u[id] = 0;
                        }

                        if(x[i] >= 0)
                        {
                                v[id] = 1;
                        }
                        else
                        {
                                v[id] = 0;
                        }
                }
        }
        
         while(t<tfinal)
	 {
	 
	 
	 t = t+k;
// initialize BCs
	if(rank_row==0)
	{
	  for(i=1;i<N_loc+1;i++){
   		 id = ID_2D(i,0,N_loc);   // bottom ghost
                u[id] = u[ID_2D(i,2,N_loc)];
                
                }
                }
	if(rank_row == q-1){
		for(i=1; i<N_loc+1; i++)
		{


                id = ID_2D(i,N_loc+1,N_loc);   // top ghost
                u[id] = u[ID_2D(i,N_loc-1,N_loc)];
                 }
                 }

	if(rank_col==0){
		for(i=1; i<N_loc+1;i++)
		{
                id = ID_2D(0,i,N_loc);   // left ghost
                u[id] = u[ID_2D(2,i,N_loc)];
                
		}
	}
	if(rank_col==q-1)
	{
		for(i=1;i<N_loc+1;i++)
		{
                id = ID_2D(N_loc+1,i,N_loc);   // right ghost
                u[id] = u[ID_2D(N_loc-1,i,N_loc)];
		
	}
}



		for(i=1; i<N_loc+1; i++)
			{
				id = ID_2D(1,i,N_loc);
				send_buffer_lft[i] = u[id];
				
			}

			for(i=1; i<N_loc+1; i++)
			{
				id = ID_2D(N_loc,i,N_loc);
				send_buffer_rgt[i] = u[id];
				
			}

			for(i=1;i<N_loc+1;i++){
     			 id = ID_2D(i,1,N_loc);
     			 send_buffer_bot[i] = u[id];
			 
    				}

  			 for(i=1;i<N_loc+1;i++){
     				 id = ID_2D(i,N_loc,N_loc);
     				 send_buffer_top[i] = u[id];
				 
    			}

			// compute neighbor indices

			 int rank_lft  = irank-1;
			 int rank_rgt  = irank+1;
			int rank_bot  = irank-q;
			int rank_top  = irank+q;

			// initialize send and recieves


    //initialize receives
    if(rank_col>0)   MPI_Irecv(recv_buffer_lft, N_loc, MPI_DOUBLE, rank_lft, 101, MPI_COMM_WORLD, &request_recv_lft);
    if(rank_col<q-1) MPI_Irecv(recv_buffer_rgt, N_loc, MPI_DOUBLE, rank_rgt, 102, MPI_COMM_WORLD, &request_recv_rgt);
    if(rank_row>0)   MPI_Irecv(recv_buffer_bot, N_loc, MPI_DOUBLE, rank_bot, 103, MPI_COMM_WORLD, &request_recv_bot);
    if(rank_row<q-1) MPI_Irecv(recv_buffer_top, N_loc, MPI_DOUBLE, rank_top, 104, MPI_COMM_WORLD, &request_recv_top);

    //initialize sends
    if(rank_col<q-1) MPI_Isend(send_buffer_rgt, N_loc, MPI_DOUBLE, rank_rgt, 101, MPI_COMM_WORLD, &request_send_rgt);
    if(rank_col>0)   MPI_Isend(send_buffer_lft, N_loc, MPI_DOUBLE, rank_lft, 102, MPI_COMM_WORLD, &request_send_lft);
    if(rank_row<q-1) MPI_Isend(send_buffer_top, N_loc, MPI_DOUBLE, rank_top, 103, MPI_COMM_WORLD, &request_send_top);
    if(rank_row>0)   MPI_Isend(send_buffer_bot, N_loc, MPI_DOUBLE, rank_bot, 104, MPI_COMM_WORLD, &request_send_bot);


//go over interior points 2:N_loc-1
			for(j=2; j<N_loc; j++){
                        for(i=2; i<N_loc; i++)
                        {

                                id = ID_2D(i,j,N_loc);
                                id_lft  = ID_2D(i-1,j,N_loc);       // index left
                                id_rgt = ID_2D(i+1,j,N_loc);        // index right
                                id_bot = ID_2D(i, j-1,N_loc);       // index bottom
                                id_top  = ID_2D(i,j+1,N_loc);       // index top

                                unew[id] = u[id] + k*f(u[id],v[id],a,b,eps) +(k/pow(h,2))*(u[id_lft]+u[id_rgt]+u[id_bot]+u[id_top]-4*u[id]);                       
                                vnew[id] = v[id] + k*g(u[id],v[id]);
                        }
                }

// enfornce no flux boundary conditions
			
			
			if(rank_row == 0){
      			for(i=1;i<N_loc+1; i++)
                        {
				j=1;
                                id = ID_2D(i,j,N_loc);
                                id_lft  = ID_2D(i-1,j,N_loc);       // index left
                                id_rgt = ID_2D(i+1,j,N_loc);        // index right
                                id_bot = ID_2D(i, j-1,N_loc);       // index bottom
                                id_top  = ID_2D(i,j+1,N_loc);       // index top

                               unew[id] = u[id] + k*f(u[id], v[id],a,b,eps) 
						+(k/pow(h,2))*(u[id_lft]+u[id_rgt]+u[id_bot]+u[id_top]-4*u[id]);
				vnew[id] = v[id] + k*g(u[id],v[id]);
                        }
                }
			

			// top boundary conditions
			if(rank_row == q-1){
			
			for(i=1; i<N_loc+1; i++)
                        {
                        j = N_loc;

                                id = ID_2D(i,j,N_loc);
                                id_lft  = ID_2D(i-1,j,N_loc);       // index left
                                id_rgt = ID_2D(i+1,j,N_loc);        // index right
                                id_bot = ID_2D(i, j-1,N_loc);       // index bottom
                                id_top  = ID_2D(i,j+1,N_loc);       // index top
				 unew[id] = u[id] + k*f(u[id], v[id],a,b,eps) 
						+(k/pow(h,2))*(u[id_lft]+u[id_rgt]+u[id_bot]+u[id_top]-4*u[id]);
				vnew[id] = v[id] + k*g(u[id],v[id]);
                        }
                }
                
                
               // left boundary

			if(rank_col == 0){
   			  for(j=1; j<N_loc+1; j++)
                        {
				i = 1;
                                id = ID_2D(i,j,N_loc);
                                id_lft  = ID_2D(i-1,j,N_loc);       // index left
                                id_rgt = ID_2D(i+1,j,N_loc);        // index right
                                id_bot = ID_2D(i, j-1,N_loc);       // index bottom
                                id_top  = ID_2D(i,j+1,N_loc);       // index top

                               unew[id] = u[id] + k*f(u[id], v[id],a,b,eps) 
						+(k/pow(h,2))*(u[id_lft]+u[id_rgt]+u[id_bot]+u[id_top]-4*u[id]);
				vnew[id] = v[id] + k*g(u[id],v[id]);
                             
                        }
                }
			
			// right boundary

			 if(rank_col==q-1){
      				for(j=1; j<N_loc+1; j++)
                        {
				
				i = N_loc;
                                id = ID_2D(i,j,N_loc);
                                id_lft  = ID_2D(i-1,j,N_loc);       // index left
                                id_rgt = ID_2D(i+1,j,N_loc);        // index right
                                id_bot = ID_2D(i, j-1,N_loc);       // index bottom
                                id_top  = ID_2D(i,j+1,N_loc);       // index top
                                 unew[id] = u[id] + k*f(u[id], v[id],a,b,eps) 
						+(k/pow(h,2))*(u[id_lft]+u[id_rgt]+u[id_bot]+u[id_top]-4*u[id]);
				vnew[id] = v[id] + k*g(u[id],v[id]);
                        }
                }
 
   //check if communication complete
    if(rank_col>0){
      MPI_Wait(&request_recv_lft,MPI_STATUS_IGNORE); //wait until left recv is complete
      for(i=1;i<N_loc+1;i++){
	id_ghost = ID_2D(0,i,N_loc); //left ghost
	u[id_ghost] = recv_buffer_lft[i]; //load receive buffer to left ghost
      }      
    }

    if(rank_col<q-1){
      MPI_Wait(&request_recv_rgt,MPI_STATUS_IGNORE); //wait until right recv is complete
      for(i=1;i<N_loc+1;i++){
	id_ghost = ID_2D(N_loc+1,i,N_loc); //right ghost
	u[id_ghost] = recv_buffer_rgt[i]; //load receive buffer to right ghost
      }      
    }

    if(rank_row>0){
      MPI_Wait(&request_recv_bot,MPI_STATUS_IGNORE); //wait until bottom recv is complete
      for(i=1;i<N_loc+1;i++){
	id_ghost = ID_2D(i,0,N_loc); //bottom ghost
	u[id_ghost] = recv_buffer_bot[i]; //load receive buffer to bottom ghost
      }      
    }
    
    if(rank_row<q-1){
      MPI_Wait(&request_recv_top,MPI_STATUS_IGNORE); //wait until top recv is complete
      for(i=1;i<N_loc+1;i++){
	id_ghost = ID_2D(i,N_loc+1,N_loc); //top ghost
	u[id_ghost] = recv_buffer_top[i]; //load receive buffer to top ghost
      }      
    }


//complete computation at the rank edges (exclude Boundaries)

		if(rank_col>0){ //compute left edge
     		 for(j=1;j<N_loc+1;j++){
		i=1;

		id = ID_2D(i,j,N_loc);
		id_lft = ID_2D(i-1,j,N_loc); //index left
		id_rgt = ID_2D(i+1,j,N_loc); //index right
		id_bot = ID_2D(i,j-1,N_loc); //index bottom
		id_top = ID_2D(i,j+1,N_loc); //index top
		
		
               unew[id] = u[id] + k*f(u[id],v[id],a,b,eps) +(k/pow(h,2))*(u[id_lft]+u[id_rgt]+u[id_bot]+u[id_top]-4*u[id]);                       
               vnew[id] = v[id] + k*g(u[id],v[id]);
      }
    }

		if(rank_col<q-1){ //compute right edge
     		 for(j=1;j<N_loc+1;j++){
		i=N_loc;

		id = ID_2D(i,j,N_loc);
		id_lft = ID_2D(i-1,j,N_loc); //index left
		id_rgt = ID_2D(i+1,j,N_loc); //index righ
		id_bot = ID_2D(i,j-1,N_loc); //index bottom
		id_top = ID_2D(i,j+1,N_loc); //index top
		
		 unew[id] = u[id] + k*f(u[id],v[id],a,b,eps) +(k/pow(h,2))*(u[id_lft]+u[id_rgt]+u[id_bot]+u[id_top]-4*u[id]);                       
                vnew[id] = v[id] + k*g(u[id],v[id]);
      }
    }
	
	if(rank_row>0){ //compute bottom edge
	      for(i=1;i<N_loc+1;i++){
		j=1;

		id = ID_2D(i,j,N_loc);
		id_lft = ID_2D(i-1,j,N_loc); //index left
		id_rgt = ID_2D(i+1,j,N_loc); //index righ
		id_bot = ID_2D(i,j-1,N_loc); //index bottom
		id_top = ID_2D(i,j+1,N_loc); //index top

		unew[id] = u[id] + k*f(u[id],v[id],a,b,eps) + (k/pow(h,2))*(u[id_lft]+u[id_rgt]+u[id_bot]+u[id_top]-4*u[id]);                       
               vnew[id] = v[id] + k*g(u[id],v[id]);
      }
    }
	if(rank_row<q-1){ //compute top edge
      	for(i=1;i<N_loc+1;i++){
		j=N_loc;
	
		id = ID_2D(i,j,N_loc);
		id_lft = ID_2D(i-1,j,N_loc); //index left
		id_rgt = ID_2D(i+1,j,N_loc); //index right
		id_bot = ID_2D(i,j-1,N_loc); //index bottom
		id_top = ID_2D(i,j+1,N_loc); //index top
		
		unew[id] = u[id] + k*f(u[id],v[id],a,b,eps) +(k/pow(h,2))*(u[id_lft]+u[id_rgt]+u[id_bot]+u[id_top]-4*u[id]);                       
                vnew[id] = v[id] + k*g(u[id],v[id]);
      }
    }
	

	// update the solution
	
	for(j=1;j<N_loc+1;j++){
      		for(i=1;i<N_loc+1;i++){
			id = ID_2D(i,j,N_loc);
			u[id] = unew[id];
			v[id] = vnew[id];
      }
    }	
}

	write_results(u, v, N_loc, nproc, irank, N, L,h);
	
	free(u);
        free(unew);
        free(v);
        free(vnew);
        free(x);
        free(y);


         return 0;
}

void getArgs(int *N, int *L, double *a, double *b, double *eps, int *tfinal, int argc, char *argv[],int irank, MPI_Comm comm)
 {
 
         if(argc != 7)
         {
                 printf("Incorect  number of commandline arguments: 7 required");
		MPI_Finalize();
		exit(1);
         }
         else
        {

                *N = atoi(argv[1]);
                *L = atoi(argv[2]);

                *a = atof(argv[3]);
                *b = atof(argv[4]);
                *eps = atof(argv[5]);
                *tfinal = atoi(argv[6]);
        }
}
	

double f(double u, double v, double a , double b, double eps)
{
        double fn =  (1/eps)*u*(1-u)*(u - (v+b)/a);
        return fn;

                        }



double g(double u, double v)
	{
	double diff = u - v;
	return diff;
	}




void write_results(double *u, double *v, int N_loc, int nproc, int irank, int N, int L, double h)
{
	int i,j, id;

  double *u_local       = malloc((N_loc)*(N_loc)*sizeof(double));
  double *u_global      = malloc((N)*(N)*sizeof(double));
  double *u_write      = malloc((N)*(N)*sizeof(double));

  double *v_local       = malloc((N_loc)*(N_loc)*sizeof(double));
  double *v_global      = malloc((N)*(N)*sizeof(double));
  double *v_write      = malloc((N)*(N)*sizeof(double));

  double x,y;



  //pack the data for gather (to avoid sending ghosts)
  int id_loc = 0;
  for(j=1;j<N_loc+1;j++){
    for(i=1;i<N_loc+1;i++){
      id = ID_2D(i,j,N_loc);
      u_local[id_loc] = u[id];
      v_local[id_loc] = v[id];
      id_loc++;
    }
  }

  MPI_Gather(u_local,id_loc,MPI_DOUBLE,u_global,id_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Gather(v_local,id_loc,MPI_DOUBLE,v_global,id_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);



  //unpack data so that it is in nice array format

  int id_write, id_global;
  int p_row, p_col;
  int q = sqrt(nproc);

  if(irank==0){
    for(int p=0; p<nproc;p++){
      p_row = p/q;
      p_col = p%q;
      for(j=0;j<N_loc;j++){
	for(i=0;i<N_loc;i++){
	  id_global = p*N_loc*N_loc + j*N_loc + i;
	  id_write  = p_row*N_loc*N_loc*q + j*N_loc*q + p_col*N_loc + i;
	  u_write[id_write] = u_global[id_global];
	  v_write[id_write] = v_global[id_global];

	}
      }
    }


     //write to file
    FILE *f;
    f = fopen("rest_barkley.dat","w"); //open file
    fprintf(f,"x, y, u, v\n");

    for(j=0; j<N; j++){
      for(i=0; i<N; i++){
	id = j*N + i;
	x = -L + i*h; //I am a bit lazy here with not gathering x and y
	y = -L + j*h;
	fprintf(f,"%f, %f, %f, %f\n",x,y,u_write[id],v_write[id]);
      }
    }
    fclose(f);
  }
  free(u_local); free(u_global); free(u_write);
  free(v_local); free(v_global); free(v_write);

}
