//=======================================================================================================

void dclear(double *a , int n);

//========================================================================================================

void dclear_matrix(double **a , int m, int n);

void dclear_matrix_int(int **a , int m, int n);
//=========================================================================================================

//===================================================================
void matvec(int nNosEl, double *p, double *ap, int neq, int nel, double ***s, int *id, int **ien);

void solucaoGmres(double ***A, double *b, double *x, int nel, int neq, double etol, int kmax, int lmax, int *id, int **ien);
