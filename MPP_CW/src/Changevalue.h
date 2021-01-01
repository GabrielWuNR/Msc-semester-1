

int isLastRow(int row, const int dims[2])  ; 
int isLastCol(int col, const int dims[2]); 
 
int typeIdx(int row, int col,  int dims[2]);
void rowcol(int rank,  int dims[2], int *row, int *col);
int  smallmap0(int row,int dims[2],int M1,int M2)  ;
int  smallmap1(int col,int dims[2],int N1,int N2) ;
 
void  random_filledsquare(int** map ,int Length,double rho);
void  SetOriginialvalue(int** old , int** smallmap,int m,int n);
void  UpdateSmallmapvalue(int** old, int** smallmap,int m,int n);
