
 /*In this source file the code realise three functions
 * intitialalltoallw: Initialize the parameters of MPI_Alltoallw
 * --sendcounts,senddispls,sendtypes,recvcounts,recvdispls,recvtypes
 *
 * mainranktoall: Initialize the parameters of MPI_Alltoallw on rank 0
 * ExchangehHalos: Parallel module for exchanging the boundary value of smallmap halo
 */
void intitialalltoallw(int size,int sendcounts[size],int senddispls[size],MPI_Datatype sendtypes[size],
                      int recvcounts[size],int recvdispls[size],MPI_Datatype recvtypes[size],
                      int smallmapsize[2]);
                      
                      
void mainranktoall(int* mapa,int M1,int M2,int N1,int N2,int mapsize[2],
                   int size,int sendcounts[size],int senddispls[size],
                   MPI_Datatype sendtypes[size],int dims[2] );
                   
                   
int ExchangehHalos(int step,int nchange,int smallmapsize[2],int ** old,int **new,
                    int ** smallmap,int left,int right ,
                    int top,int down,MPI_Comm topo_comm,int rank );





