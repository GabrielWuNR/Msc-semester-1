#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "src_h/percolate.h"
#include "src_h/arralloc.h"
#include "src_h/Changevalue.h"


 
void intitialalltoallw(int size,int sendcounts[size],int senddispls[size],MPI_Datatype sendtypes[size],
                      int recvcounts[size],int recvdispls[size],MPI_Datatype recvtypes[size],
                      int smallmapsize[2])
{ int i;
  for (  i=0; i<size; i++)
     {
        recvcounts[i] = 0;
        recvdispls[i] = 0;
        recvtypes[i] = MPI_INT;

        sendcounts[i] = 0;
        senddispls[i] = 0;
        sendtypes[i] = MPI_INT;
      }
    
    recvcounts[0] = smallmapsize[0]*smallmapsize[1];
    recvdispls[0] = 0;
    

}

void mainranktoall(int* mapa,int M1,int M2,int N1,int N2,int mapsize[2],int size,int sendcounts[size],int senddispls[size],MPI_Datatype sendtypes[size],int dims[2])
{ int i,j;
 MPI_Datatype blocktypes[4];
    int subsizes[2];
    int starts[2] = {0,0};
	 for (i = 0; i < 2; i++) 
    {
		  if(i%2==0)subsizes[0] =M1;
		  if (i % 2 !=0)subsizes[0] = M2;
		for (j = 0; j < 2; j++)
     {
			if (j % 2 == 0)subsizes[1] = N1;
			if (j % 2 != 0)subsizes[1] = N2; 
		 
     MPI_Type_vector(subsizes[0],subsizes[1], mapsize[1], MPI_INT,&blocktypes[2 * i + j]);//Create 2D Vector Derived Data Type
     MPI_Type_create_resized(blocktypes[2 * i + j], 0, subsizes[0]*subsizes[1]*sizeof(int), &blocktypes[2 * i + j]);
	 	 MPI_Type_commit(&blocktypes[2 * i + j]);
	  	}
    }  
  /*
   * According to the characteristics of the Cartesian coordinate system, 
   * the code uses the for loop to judge the progress of each level 
   * when it is in the last row or the last column.
   */
     for (  i=0; i<size; i++)
        {
          int row, col; 
          int smallmapsize0,smallmapsize1;
          rowcol(i, dims, &row, &col);
          smallmapsize0= smallmap0(row,dims,M1,M2) ;
          smallmapsize1= smallmap1(col,dims,N1,N2) ;
          sendcounts[i] = 1;
          senddispls[i] = (row*M1*mapsize[1] + col*N1)*sizeof(int);
          int idx = typeIdx(row, col, dims);
          sendtypes[i] = blocktypes[idx];     
        }

}

int ExchangehHalos(int step,int nchange,int smallmapsize[2],int ** old,int **new,int ** smallmap,int left,int right ,int top,int down,MPI_Comm topo_comm,int rank )
{  MPI_Status status1,status2,status3,status4;
   MPI_Request request1,request2,request3,request4;
// maxstep = 16*L;
step = 1;
 int printfreq = 100;
int maxstep=16*L;
 
  int i,j;
  
  
  
  while (step <= maxstep)
    { 
     nchange = 0;  
   for (j=1; j <= smallmapsize[1]; j++)
  {
    old[0][j]   = old[smallmapsize[0]][j];
    old[smallmapsize[0]+1][j] = old[1][j];
  } 
  
    MPI_Issend(&old[smallmapsize[0]][1],smallmapsize[1],MPI_INT,right,0,MPI_COMM_WORLD,&request3);//send the right halo to right neighbour
          MPI_Recv(&old[0][1],smallmapsize[1],MPI_INT,left,0,MPI_COMM_WORLD,&status3);//receive the left halo from left neighbour
            MPI_Wait(&request3, &status3);
       MPI_Issend(&old[1][1],smallmapsize[1],MPI_INT,left,0,MPI_COMM_WORLD,&request4);//send the right halo to right neighbour
        MPI_Recv(&old[smallmapsize[0]+1][1],smallmapsize[1],MPI_INT,right,0,MPI_COMM_WORLD,&status4);//receive the right halo from the right neighbour      
         MPI_Wait(&request4, &status4);
   
    for(i=1;i<=smallmapsize[0];i++)
      {
         MPI_Issend(&old[i][smallmapsize[1]],1,MPI_INT,top,0,topo_comm,&request1);//send the right halo to right neighbour 
          MPI_Recv(&old[i][0],1,MPI_INT,down,0,topo_comm,&status1);//receive the left halo from left neighbour
           MPI_Wait(&request1, &status1); 
         MPI_Issend(&old[i][1],1,MPI_INT,down,0,topo_comm,&request2);//send the left halo to left neighbour
             MPI_Recv(&old[i][smallmapsize[1]+1],1,MPI_INT,top,0,topo_comm,&status2);//receive the right halo from the right neighbour
            MPI_Wait(&request2, &status2);            
      }  
  
      int oldval,newval;
   
    for (i=1; i<=smallmapsize[0]; i++)
	{
	  for (j=1; j<=smallmapsize[1]; j++)
	    {
	      oldval = old[i][j];
	      newval = oldval;
 
	      if (oldval != 0)
		{
		  if (old[i-1][j] > newval) newval = old[i-1][j];
		  if (old[i+1][j] > newval) newval = old[i+1][j];
		  if (old[i][j-1] > newval) newval = old[i][j-1];
		  if (old[i][j+1] > newval) newval = old[i][j+1];

		  if (newval != oldval)
		    {
		      ++nchange;
		    }
		}
	      new[i][j] = newval;
	    }
	}
 
    int judgestop,judgestop1;
   MPI_Allreduce(&nchange,&judgestop,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);  
   if(judgestop==0)   break;
   
   
   
       int iter;
 
      
     int mapall,average=0;
    for(i=1;i<=smallmapsize[0];i++)
    {
       for(j=1;j<=smallmapsize[1] ;j++)
     {
      average+=old[i][j];
      }
    }
    
    
    int map_ave;
   MPI_Allreduce(&average,&mapall,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);  
    
     map_ave =mapall/(L*L);
     
               
    if (step % printfreq == 0)
     {  
    
            
        if(rank==0)
        {  
          printf("percolate: number of changes on step %d is %d, rank is %d and average number of array is %d \n",
         step, nchange,rank,map_ave);
         }
         
     }  
     
     
      for (i=1; i<=smallmapsize[0]; i++)
	    {
	  for (j=1; j<=smallmapsize[1]; j++)
	      {  
        old[i][j] = new[i][j];
         }
    }
      step++;
      
    }
    if(rank==0)
    {
     printf("\nstep is %d!\n",step);
    }
    return nchange;  
   }
  
  
  //void sendBorder()  