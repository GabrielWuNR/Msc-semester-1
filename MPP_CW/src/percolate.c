#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "arralloc.c"
#include "src_h/arralloc.h"
#include "src_h/percolate.h"
#include "src_h/Changevalue.h"
#include "src_h/parallel.h"
 

int main(int argc, char *argv[])


{
 
  /*
   *  Define the main arrays for the simulation
   */
    int **old , **new ;
    
    
  /*
   *  Additional array WITHOUT halos for initialisation and IO. This
   *  is of size LxL because, even in our parallel program, we do
   *  these two steps in serial
   */
    int **map ;
    int **smallmap ;
    
  /*
   *  Variables that define the simulation
   */
    int seed;
    double rho; 
  /*
   *  Local variables
   */
   //  int i, j, step, maxstep, oldval, newval, nchange, printfreq;
     int i, j, step, nchange;
  
   
  /*
   *  M1,M2,N1,N2: Variables that define the width and length of smallmap and map 
   *  M1:When L cannot divide the number of processes per row, the value of M is rounded down. Same as N1.
   *  M2: M value in the last row   N2:N value in the last column
   *  per_row_size,per_col_size:variables that define how many processes in one row or in one col
   *   per_row_size:  Number of processes per row in the topology
   *  per_col_size : //Number of processes per column in the topology 
   */ 
   
   
    int M1,N1,M2,N2;
    int per_row_size;
    int per_col_size;
    int mapsize[2];
       
    int smallmapsize[2];
  
   mapsize[0]=L;
   mapsize[1]=L;
  
   //initialize MPI
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);  //Get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  //Get the rank of the process

   /*Create a Cartesian coordinate system 
    *based on the original communication domain
    */
     
   int nnodes;
   int ndims;//Topological dimension
   int dims[2] = {0,0};
    int peroid[2]={1,0};
   int new_rank;
   int coords[2]={0,0};
   int left,right;
   int top,down;
    MPI_Comm  topo_comm;  
  
   nnodes=size;
   ndims=2;
 
   
  
   MPI_Dims_create(size,2,dims); 
   MPI_Cart_create(MPI_COMM_WORLD,2,dims,peroid,0,&topo_comm);
   MPI_Comm_rank(topo_comm, &new_rank);
   MPI_Cart_coords(topo_comm,new_rank, 2, coords);  
   MPI_Cart_shift(topo_comm, 0, 1,&left,&right);	// left& right neighbours	  
   MPI_Cart_shift(topo_comm, 1, 1,&down,&top);	//top & down neighbours*/
  
  
   
   
   per_row_size=dims[0];
   per_col_size=size/per_row_size;
   M1=(int) L/ per_row_size;
   N1=(int) L/ per_col_size;
   M2=L-M1*(per_row_size-1);
   N2=L-N1*(per_col_size-1);
   
   
   
   //M of smallmap 
   smallmapsize[0] =M1;
    //N of smallmap
   smallmapsize[1] = N1; 
   /*Determine if the process is in the last row in the topology. 
    *If so, update its smallmap M value to M2*/
   if ( isLastRow(coords[0],dims)==1)  smallmapsize[0]=M2;
   /*Determine if the process is in the last column in the topology.
   * If so, update its smallmap N value to N2*/
   if ( isLastCol(coords[1],dims)==1) smallmapsize[1]=N2;
     
   // smallmapsize(dims,coords,M1,M2, N1,N2,smallmapsize);
  
    
    
    
   //initialise smallmap
   smallmap = arralloc(sizeof(int), 2, smallmapsize[0], smallmapsize[1]); 
    //initialise  map
    map = arralloc(sizeof(int), 2, L, L); 
     //initialise smallmap with halo
    old = arralloc(sizeof(int), 2, smallmapsize[0]+2 ,smallmapsize[1]+2);  
    new = arralloc(sizeof(int), 2, smallmapsize[0]+2 ,smallmapsize[1]+2);  
    
    
    
   int *mapa=NULL;//map start address
     
 
 
 
   if (argc != 2)
    {
      printf("Usage: percolate <seed>\n");
      return 1;
    }
  
   rho = 0.411;
   seed = atoi(argv[1]);

   rinit(seed);
 
  
 
  if(size != P)
  {
  printf("the size is not equal to P, program break!");
  MPI_Finalize();
  exit(0);
  }
  
  int times;
  double time_aver[3],time_aver_step;
  
  /*
   *  Initialise map with density rho. Zero indicates rock, a positive
   *  value indicates a hole. For the algorithm to work, all the holes
   *  must be initialised with a unique integer
   */

   if(rank == 0)
   {
   random_filledsquare(map,L,rho);
   }
  
  
  
   
    
    /*
    *Send map data from the main process via MPI_Alltoallw
    */
    int starts[2]   = {0,0};               
    int sendcounts[size];//Number of messages sent to each process
    int senddispls[size];//Message displacement sent to each process
    MPI_Datatype sendtypes[size];//Message data type sent to each process
    MPI_Request request;
    MPI_Status status;
 
    
    int recvcounts[size];//Number of messages received by each process
    int recvdispls[size];//Message displacement received by each process
    MPI_Datatype recvtypes[size];//Message data type received by each process
    
   /*We only need to specify the data type and data size of
    *each process to be sent in the main process.
    *The other processes only need to set the send message count to 0.
    */    
    intitialalltoallw(size,sendcounts,senddispls,sendtypes,
    recvcounts,recvdispls,recvtypes,smallmapsize);
   
    
   /* According to the case where both sides of the map cannot be divisible,
    * there are up to four MP_Vectors of different sizes,
    * the sizes are M1 * N1, M1 * N2, M2 * N1, and M2 * N2.
    */
   
      if (rank == 0)
  {
    mapa = &(map[0][0]);
    
    mainranktoall(mapa,M1,M2,N1,N2,mapsize,size,
    sendcounts,senddispls,sendtypes,dims);
    
    
  }
       
  //calculate the sending time of smallmap
    MPI_Barrier(topo_comm);
   double start0 = MPI_Wtime();
    MPI_Alltoallw(mapa, sendcounts, senddispls, sendtypes,
              &(smallmap[0][0]), recvcounts, recvdispls, recvtypes, 
              topo_comm); 
 
      double end0 = MPI_Wtime();
     double time0= end0-start0;   
     time_aver[0]+=time0;   
  MPI_Barrier(topo_comm);
  
  
    
  /*
   * Initialise the old array: copy the LxL array map to the centre of
   * old, and set the halo values to zero.
   */ 
   SetOriginialvalue(old,smallmap,smallmapsize[0],smallmapsize[1]);

  
  //step = 1;
  nchange = 1;
  
  
  MPI_Barrier(topo_comm);  
  double start = MPI_Wtime();


/*
  Using the while loop to exchange the boundary value of the smallmap 
  with the halo to the upper, lower, left, and right processes, 
  and update the value of smallmap itself
  */
  
  nchange=ExchangehHalos(step,nchange,smallmapsize,old,new,smallmap,left,right ,top,down,topo_comm,rank );
 
    
   double end = MPI_Wtime();
  
   time_aver[1]+=(end-start);
    
   
    if(rank==0)
    {
     printf("\n\n The time cost in while is %f!\n\n\n",end-start);
    }
   MPI_Barrier(topo_comm);
    
    
    
  if (nchange != 0)
    {
     printf("percolate: WARNING max steps = %d reached before nchange = 0\n",
	    nchange);
    }
    
    
  UpdateSmallmapvalue( old , smallmap,smallmapsize[0],smallmapsize[1]);
  
  
  //calculate the time of receiving the smallmap to map in main rank 0
    MPI_Barrier(topo_comm);
    
    double start1 = MPI_Wtime();
    
    MPI_Alltoallw( &(smallmap[0][0]), recvcounts, recvdispls, recvtypes,
                 mapa, sendcounts, senddispls, sendtypes, 
                  topo_comm); 
                  
     double end1 = MPI_Wtime();
     double time1= end1-start1;   
     time_aver[2]+=time1;           
     MPI_Barrier(topo_comm);
   
    
    
 free(smallmap);
 
  //Judge whether map has percolated in main rank 0
  if(rank == 0)
  {
    JudgePercolate(map);
  percwritedynamic("map.pgm", map,L, 9);
 }
 
  
 //print the time of three parallel mode 
   if(rank==0){
    printf("\n\n\n\n average time of sending smallmap is %f  \n",time_aver[0]);
    printf("average time in switching the halo is %f\n",time_aver[1]);
    printf("average time of receiving smallmap is %f  \n",time_aver[2]);
   } 
     MPI_Finalize();
 
  return 0;
}


