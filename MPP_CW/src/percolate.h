 
/*
 *  Main header file for percolation code.
 */
/*
 * Factor for number of processors
 */
 
/* 
 *  System size L
 */
 
#define L 1000
 
/*
 *  Although overall system is square, i.e. size L x L, we will define
 *  different variables for the first and second dimensions. This is
 *  because, in the parallel code, the local arrays will not be
 *  square. For example, using a simple 1D decomposition over P
 *  processes, then M = L/P and N = L
 */
#define P 15

 
/*
 *  Prototypes for supplied functions
 */

/*
 *  Visualisation
 */
 



 
void percwritedynamic(char *percfile, int **map, int l, int ncluster);
/*
 *  Random numbers
 */
int * dividenum(int num);
void rinit(int ijkl);
float uni(void);
 
void JudgePercolate(int** map );
 
 

 
 
