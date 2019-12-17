
#include<stdio.h>
#include<stdlib.h>

//===================================================================
//funcao que aloca dinamicamente uma vetor de doubles
double *criaVetorD(int n){
	double *v;   //ponteiro para vetor

	if (n < 1){
		printf ("\n** Erro: Parametro invalido **\n");
		return NULL;
	}
	//aloca o vetor
	v = (double*)malloc(n*sizeof(double));// Um vetor de n doubles
	if (v == NULL){
		printf ("** Erro: Memoria Insuficiente **");
		return NULL;
	}
	return v;
}

//===================================================================
//funcao que aloca dinamicamente uma vetor de inteiros
int *criaVetorI(int n){
	int *v;   //ponteiro para vetor

	if (n < 1){
		printf ("\n** Erro: Parametro invalido **\n");
		return NULL;
	}
	//aloca o vetor
	v = (int*)malloc(n*sizeof(int));// Um vetor de n inteiros
	if (v == NULL){
		printf ("** Erro: Memoria Insuficiente **");
		return NULL;
	}
	return v;
}

//===================================================================
//funcao que aloca dinamicamente uma matriz de doubles
double **criaMatrizD(int m, int n){    /* float mxn-matriz */
	int i;
	double **a;

	a = (double **) calloc(m,sizeof(double *));
	if (a == NULL){
		printf ("** Erro: Memoria Insuficiente **");
		return NULL;
		//fprintf(stderr,"Nao pode alocar memoria");
		//exit(1);
	}
	for(i = 1; i <= m; i++){
		a[i-1] = (double *) calloc(n, sizeof(double));
		if(a[i-1] == NULL){
			printf ("** Erro: Memoria Insuficiente **");
			return NULL;
			//fprintf(stderr,"Nao pode alocar memoria");
			//exit(1);
		}
	}
	return a;
}

//===================================================================
//funcao que aloca dinamicamente uma matriz de inteiros
int **criaMatrizI(int m, int n){    /* int mxn-matriz */
	int i;
	int **a;

	a = (int **) calloc(m,sizeof(int *));
	if (a == NULL){
		printf ("** Erro: Memoria Insuficiente **");
		return NULL;
		//fprintf(stderr,"Nao pode alocar memoria");
		//exit(1);
	}
	for(i = 1; i <= m; i++){
		a[i-1] = (int *) calloc(n, sizeof(int));
		if(a[i-1] == NULL){
			printf ("** Erro: Memoria Insuficiente **");
			return NULL;
			//fprintf(stderr,"Nao pode alocar memoria");
			//exit(1);
		}
	}
	return a;
}

//===================================================================
//funcao que libera espaco reservado para a matriz de double
double **destruirMatrizD(double **mat, int numlinhas, int numcolunas){
	int  i;

	if (mat == NULL)
		return NULL;
	for (i=0; i<numlinhas; i++)
		free(mat[i]);    //libera as linhas da matriz
	free(mat);         //libera a matriz (vetor de ponteiros)

	return NULL;        //retorna um ponteiro nulo
}

//===================================================================
//funcao que libera espaco reservado para a matriz de int
int **destruirMatrizI(int **mat, int numlinhas, int numcolunas){
	int  i;
	if (mat == NULL)
		return NULL;
	for (i=0; i<numlinhas; i++)
		free(mat[i]);    //libera as linhas da matriz
	free(mat);         //libera a matriz (vetor de ponteiros)
	return NULL;        //retorna um ponteiro nulo
}

//===================================================================

void imprimirVetorD(double *x, int n){
	int i;
	printf("\nvetor = [");
	for(i = 0; i < n; i++){
		printf("%.4lf",x[i]);
		if(i == n-1)
			printf("]");
		else
			printf(", ");
	}
	printf("\n");
}

//===================================================================

void imprimirVetorI(int *x, int n){
	int i;
	printf("\nvetor = [");
	for(i = 0; i < n; i++){
		printf("%d, ",x[i]);
		if(i == n-1)
			printf("]");
		else
			printf(", ");
	}
	printf("\n");
}
//===================================================================

void  imprimirMatrizD(double **mat, int m, int n){
	int i, j;
	//printf("\nm := matrix(12,12,[");
	printf("\n\n");
	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			printf("%.2lf,   ",mat[i][j]);
			if(j == n-1)
				printf("\n");
		}
	}
	//printf("]);\n");
}

//===================================================================

void  imprimirMatrizI(int **mat, int m, int n){
	int i, j;
	//printf("\nm := matrix(12,12,[");
	printf("\n\n");
	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			printf("%d,   ",mat[i][j]);
			if(j == n-1)
				printf("\n");
		}
	}
	//printf("]);\n");
}
//=====================================================================
double ***criaMatrizTD(int p, int m, int n){
	int i, j;
	double ***a;

	a = (double ***) calloc(p,sizeof(double **));
	if (a == NULL){
		printf ("** Erro: Memoria Insuficiente **");
		return NULL;
	}
	for (i=0; i<p; i++){
		a[i] = (double **) calloc(m,sizeof(double *));
		for (j=0; j<m; j++){
			a[i][j] = (double *) calloc(n,sizeof(double));
		}
	}
	return a;
}
