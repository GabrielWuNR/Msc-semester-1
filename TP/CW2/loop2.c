/* B157092 */
#include <stdio.h>
#include <math.h>


#define N 729
#define reps 1000 
#include <omp.h> 

double a[N][N], b[N][N], c[N];
int jmax[N];  


void init1(void);
void init2(void);
void runloop(int); 
void loop1chunk(int, int);
void loop2chunk(int, int);
void valid1(void);
void valid2(void);


  /* Affinity schedule functions */
  
  
  
  
 
int DetectMostLoadThread(int ,int ,int* ,int* );
  
void UpdateStop_Start(int ,int* ,int* ,int ,int *,int *);
 
int main(int argc, char *argv[]) { 

  double start1,start2,end1,end2;
  int r;

  init1(); 

  start1 = omp_get_wtime(); 

  for (r=0; r<reps; r++){ 
    runloop(1);
  } 

  end1  = omp_get_wtime();  

  valid1(); 

  printf("Total time for %d reps of loop 1 = %f\n",reps, (float)(end1-start1)); 


  init2(); 

  start2 = omp_get_wtime(); 

  for (r=0; r<reps; r++){ 
    runloop(2);
  } 

  end2  = omp_get_wtime(); 

  valid2(); 

  printf("Total time for %d reps of loop 2 = %f\n",reps, (float)(end2-start2)); 

} 

void init1(void){
  int i,j; 

  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      a[i][j] = 0.0; 
      b[i][j] = 3.142*(i+j); 
    }
  }

}

void init2(void){ 
  int i,j, expr; 

  for (i=0; i<N; i++){ 
    expr =  i%( 3*(i/30) + 1); 
    if ( expr == 0) { 
      jmax[i] = N;
    }
    else {
      jmax[i] = 1; 
    }
    c[i] = 0.0;
  }

  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      b[i][j] = (double) (i*j+1) / (double) (N*N); 
    }
  }
 
} 


void runloop(int loopid) 
 { 

   int nthreads = omp_get_max_threads();  
    int ipt = (int) ceil((double)N/(double)nthreads);  
    float per_fraction=1.0/nthreads;
    int iter_start[nthreads];
     int iter_stop[nthreads];
   
  #pragma omp parallel default(none) shared(loopid,nthreads,ipt,iter_start,iter_stop) 
  {
  
   int myid  = omp_get_thread_num(); 
    int lo = myid*ipt; 
    int hi = (myid+1)*ipt; 
    if (hi > N) hi = N;  
    iter_start[myid]=lo;
    iter_stop[myid]=hi;
    
    int re_excuteJ=0;
    int iStart=0, iStop=0;
    
    #pragma omp barrier 
    while(re_excuteJ!=-1)
    {
      switch (loopid) 
       { 
       case 1: loop1chunk(iStart,iStop); break;
       case 2: loop2chunk(iStart,iStop); break;
       }
    /*Judge whether all the processes has completed its own work,
    *if not, Update its Upper and lower bounds
    
    */
    #pragma omp critical
			{
      re_excuteJ=DetectMostLoadThread(nthreads,myid,iter_start,iter_stop);
      
       UpdateStop_Start(re_excuteJ,iter_start,iter_stop,nthreads,&iStart,&iStop);
      }
     
    } 
  }
}
 /*Determine whether the thread is idle. 
   *If idle, the function returns the thread
   * which needs to be shared 1/p of 
   *the heaviest workload in all threads.
   */
int DetectMostLoadThread(int nthread,int threadId,int* iter_start,int* iter_stop)
{
  int current_position;
   if(iter_stop[threadId]-iter_start[threadId]>0)
   {
     current_position=threadId;
   }
   else
   {
   /*Used to determine the currently loaded thread number, initialized to -1
    *If all processes have been executed, 
    *the final return value is MostLposition initial value -1
    */
     int MostLposition=-1;
     int gap_initial=0;
     int i;
     for(i=0;i<nthread;i++)
     {
       int gap=iter_stop[i]-iter_start[i];
       if(gap>gap_initial)
       {
         MostLposition=i;
         gap_initial=gap;
       }
     }
     current_position=MostLposition;        
   }
   
  return current_position;
}
 /*Update the workload of the for loop 
    *shared by the current thread, which means updating
    *the upper and lower bounds it undertakes in the for loop
    */

 void UpdateStop_Start(int MostLposition,int* iter_start,int* iter_stop,int nthread,int *iStart,int *iStop)
{ 
  if(MostLposition!=-1)
  {
   float per_fraction=1.0/nthread;
   
   int remain_iter=iter_stop[MostLposition]-iter_start[MostLposition];
   int newgap = (int)ceil(per_fraction * remain_iter);
  
  	*iStart = iter_start[MostLposition];
	  *iStop = iter_start[MostLposition] + newgap;
	  iter_start[MostLposition] =iter_start[MostLposition] + newgap;
   }
 
}




void loop1chunk(int lo, int hi) { 
  int i,j; 
  
  for (i=lo; i<hi; i++){ 
    for (j=N-1; j>i; j--){
      a[i][j] += cos(b[i][j]);
    } 
  }

} 



void loop2chunk(int lo, int hi) {
  int i,j,k; 
  double rN2; 

  rN2 = 1.0 / (double) (N*N);  

  for (i=lo; i<hi; i++){ 
    for (j=0; j < jmax[i]; j++){
      for (k=0; k<j; k++){ 
	c[i] += (k+1) * log (b[i][j]) * rN2;
      } 
    }
  }

}

void valid1(void) { 
  int i,j; 
  double suma; 
  
  suma= 0.0; 
  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      suma += a[i][j];
    }
  }
  printf("Loop 1 check: Sum of a is %lf\n", suma);

} 


void valid2(void) { 
  int i; 
  double sumc; 
  
  sumc= 0.0; 
  for (i=0; i<N; i++){ 
    sumc += c[i];
  }
  printf("Loop 2 check: Sum of c is %f\n", sumc);
} 
