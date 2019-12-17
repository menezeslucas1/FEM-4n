//===================================================================
//funcao que aloca dinamicamente uma vetor de doubles
double *criaVetorD(int n);

//===================================================================
//funcao que aloca dinamicamente uma vetor de inteiros
int *criaVetorI(int n);

//===================================================================
//funcao que aloca dinamicamente uma matriz de doubles
double **criaMatrizD(int m, int n);

//===================================================================
//funcao que aloca dinamicamente uma matriz de inteiros
int **criaMatrizI(int m, int n);
//===================================================================
//funcao que libera espaco reservado para a matriz de double
double **destruirMatrizD(double **mat, int numlinhas, int numcolunas);

//===================================================================
//funcao que libera espaco reservado para a matriz de int
int **destruirMatrizI(int **mat, int numlinhas, int numcolunas);

//===================================================================

void  imprimirVetorD(double *x, int n);

//===================================================================

void imprimirVetorI(int *x, int n);

//===================================================================

void  imprimirMatrizD(double **mat, int m, int n);

//===================================================================

void  imprimirMatrizI(int **mat, int m, int n);

double ***criaMatrizTD(int p, int m, int n);
