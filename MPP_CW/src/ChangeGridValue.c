#include <stdio.h>
#include <stdlib.h>
#include "src_h/percolate.h"

 

 
 /*
   *  Initialise map with density rho. Zero indicates rock, a positive
   *  value indicates a hole. For the algorithm to work, all the holes
   *  must be initialised with a unique integer
   */
void  random_filledsquare(int** map ,int Length,double rho)
{
  int nhole = 0;
  int i,j;
   double r;

  for (i=0; i < Length; i++)
    {
      for (j=0; j < Length; j++)
	{
	  r=uni();
	  
	  if(r < rho)
	    {
	      map[i][j] = 0;
	    }
	  else
	    {
	      nhole++;
	      map[i][j] = nhole;
	    }
	}
    }
 
 printf("percolate: rho = %f, actual density = %f\n",
 rho, 1.0 - ((double) nhole)/((double) Length*Length) );
}


 /*
   * Initialise the old array: copy the LxL array map to the centre of
   * old, and set the halo values to zero.
   */
  
void  SetOriginialvalue(int** old , int** smallmap,int m,int n)
{ int i,j;
for (i=1; i <= m; i++)
    {
      for (j=1; j <= n; j++)
	{
	  old[i][j] = smallmap[i-1][j-1];
	}
    }
   
   for (i=0; i <= m+1; i++)  // zero the bottom and top halos
    {
      old[i][0]   = 0;
      old[i][n+1] = 0;
    }
 
}

/*
   * Update smallmap value using the value of old array
   * 
   */
void  UpdateSmallmapvalue(int** old, int** smallmap,int m,int n)
{ int i,j;
for (i=1; i<=m; i++)
    {
      for (j=1; j<=n; j++)
	{
	  smallmap[i-1][j-1] = old[i][j];
	}
    }

}


/*
   * Judge whether the process located 
   *in the last row of Cartesian coordinate system 
   */


 int isLastRow(int row, const int dims[2]) 
{   if(row== dims[0]-1)
   { return 1;}
    else
    { return 0;}
}
/*
   * Judge whether the process located 
   *in the last column of Cartesian coordinate system 
   */

int isLastCol(int col, const int dims[2])
{   if(col== dims[1]-1)
   { return 1;}
    else
    { return 0;}
}
 /*
   * Judge the blocktype of smallmap in each processes
   */

int typeIdx(int row, int col,  int dims[2]) {
    int lastrow = (row == dims[0]-1);
    int lastcol = (col == dims[1]-1);

    return lastrow*2 + lastcol;
}
  /*
   * On the main process, determine 
   *which row the process with level I is in
   */

void rowcol(int rank,  int dims[2], int *row, int *col) {
    *row = rank/dims[1];
    *col = rank % dims[1];
}

   /*
   * On the main process, determine 
   *the value of M in  smallmap on each process
   */
    
int  smallmap0(int row,int dims[2],int M1,int M2) 
{ 
  int width=M1;
   int islastRow;
   islastRow=isLastRow(row,dims);
    if ( islastRow==1) 
 {  width=M2;
   
 }
 return width;
}

  /*
   * On the main process, determine 
   *the value of N in  smallmap on each process
   */
int  smallmap1(int col,int dims[2],int N1,int N2) 
{
 int length=N1;
   int islastCol;
   islastCol=isLastCol(col,dims);
    if ( islastCol==1) 
 {  length=N2;
   
 }
 return length;
}

	      

 
   
