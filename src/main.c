#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include "matriz.h"
#include "malha.h"
#include "gmres.h"
#include "lcd.h"
#include "graphics.h"
#include "parametro.h"
#include "metgalerkin_CN.h"

int ITER_GMRES;
int ITER_LCD;

void gravar(double **COORD,//matriz de coordenadas
			double *solG,
			int nnos,//número de nós na malha
			int nosx,//número de nós no eixo x
			int nosy,// e no eixo y(malha retangular)
			int cont_t
			)
{
	FILE *oout;
	char str[100];

	mkdir("resultados", 0755);
	mkdir("resultados/dados", 0755);
	sprintf(str, "resultados/dados/resultado%d.txt",cont_t);
	oout = fopen(str,"w");
	int j=0;
	for(int i = 0; i < nnos; i++){
		fprintf(oout,"%2.13lf  %2.13lf  %lf\n",COORD[i][0],COORD[i][1],solG[i]);
		j++;
		if (j == nosx)
		{
			fprintf(oout,"\n");
			j=0;
		}
	}
	fclose(oout);
}

int main(){
	int i, j;
	int nosx, nosy;//número de nós no eixo x e no eixo y(malha retangular)
	int nnos;//número de nós na malha
	int nel;//número de elementos na malha
	int neq;//número de equações
	int cc;//tipo da fronteira (com nó prescrito ou não)
	int tipoMalha;//tipo da malha
	int nNosEl;//número de nós do elemento
	double limInfx, limSupx, limInfy, limSupy;//limites superiores e inferiores para x e para y
	double delta_t;//variação no tempo
	double tempo;//instante de tempo atual
	double tfim;//tempo final
	double alpha;
	double *sol;//solução do sistema linear
	double *sol_ant;//solução no tempo anterior
	double *solG;
	double *BV;//condição de fronteira (BOUNDCOUND)
	double **COORD;//matriz de coordenadas
	double ***A;//matriz do sistema linear A*sol = P
	double ***AA;//(M - delta_t*K/2)
	double *F;//Vetor de força
	double *F_ant;
	double *P;//vetor do sistema linear A*sol = P
	double tol;
	int *ID;//identifica a equação associada a cada nó não prescrito
//	int NLF[nel][2];//informa as arestas que o elemento possui na fronteira
	int *ND; //indica fronteira de Neumann
	int **HN; //indica as arestas de Neumann
	int **IEN;//matriz de conectividade
	int kmax, lmax, cont_t;
	clock_t tInicio, tFinal;//horrio de início e término do cálculo
	double tGasto;//tempo total gasto calculando
	FILE *iin;//arquivo de entrada

	//******************************************************
	// ... Entrada de Dados
	//******************************************************
	tInicio = clock();
	iin = fopen("input.txt", "r");
	char ch;
	while ((ch=fgetc(iin))!='\n');//ignora a primeira linha do arquivo de entrada
	fscanf(iin,"%d  %d  %d  %lf  %lf  %lf  %lf %d\n",&nosx,&nosy,&tipoMalha,&limInfx,&limSupx,&limInfy,&limSupy,&cc);

//	nel = ((nosx - 1)*(nosy - 1))*2;
	nel = ((nosx - 1)*(nosy - 1));
	printf("nel=%d\n",nel);
	nnos = nosx*nosy;
//	Aresta Front[nnos];//indica pontos na fronteira
//	int **NLF;//informa as arestas que o elemento possui na fronteira

	//******************************************************
	// ... Tipo de Fronteira
	//******************************************************


	//******************************************************
	// ... Malha
	//******************************************************
	nNosEl = 4;
	printf("\n\n ... Gerando a malha\n");
//	IEN = criaMatrizI(nel,3);
	IEN = criaMatrizI(nel,nNosEl);
	COORD = criaMatrizD(nnos,2);
	ID = criaVetorI(nnos); //identifica a equacao associada a cada noh nao prescrito
	ND = criaVetorI(nnos);
	HN = criaMatrizI(nnos,2);
//	Front = criaVetorI(nnos);
//	NLF = criaMatrizI(nel,2);
	neq = condFront(ID,ND,HN,nnos,nosx,nosy,cc);
	gerarMalha(nosx,nosy,tipoMalha,limInfx,limSupx,limInfy,limSupy,IEN,COORD);
//printf("neq=%d depois\n",neq);

/*
	for(int asdf=0;asdf<nosx;asdf++){
		printf("\n");
		for(int qwer=0;qwer<nosy;qwer++)
			printf("%3d ",ID[asdf*nosx+qwer]);
	}
*/
	BV = criaVetorD(nnos);
	lerCondFront(BV,ID,nosx,nosy);
/*
	for (int asdf=0; asdf<nnos; asdf++)
		printf("%.3f ", BV[asdf]);
*/
	A = criaMatrizTD(nel,nNosEl,nNosEl);
	AA = criaMatrizTD(nel,nNosEl,nNosEl);
	F = criaVetorD(neq);
	F_ant = criaVetorD(neq);
	P = criaVetorD(neq);

	sol = criaVetorD(neq);
	sol_ant = criaVetorD(neq);

	//solucao inicial
	solucaoInicial(sol,ID,COORD,nnos);
	for(i = 0; i < neq; i++)
    {
        F[i]=0;
    }

///////////////
	solG = criaVetorD(nnos);
	j = 0;
	for(i = 0; i < nnos; i++){
		if(ID[i] != 0){
			solG[i] = sol[j];
			j++;
		}
		else
			solG[i] = BV[i];
	}
	gravar(COORD, solG, nnos, nosx, nosy, 0);
/*
	j=0;
	oout = fopen("galerkin0.txt","w");
	for(i = 0; i < nnos; i++){
		fprintf(oout,"%2.13lf  %2.13lf  %lf\n",COORD[i][0],COORD[i][1],solG[i]);
		if (j==nosx-1)//numero de nos
		{
			fprintf(oout,"\n");
			j=0;
		}
		else j++;
	}
	fclose(oout);
*/
	free(solG);

//	scriptFigGnuplot(nosx,nosy,limInfx,limSupx,limInfy,limSupy,0);
//////////////
	kmax = kmax_gmres();
	lmax = lmax_gmres();
	tol = tol_gmres();
	ITER_GMRES = 0;
	ITER_LCD = 0;

	cont_t = 0;
	delta_t = 0.0025;
	tempo = 0.0;
	tfim = 200*delta_t;
	alpha = 0.5;
	//PROCESSO ITERATIVO NO  TEMPO - ALGORITMO CRANK-NICOLSON
	do{
		cont_t = cont_t + 1; //contador de tempo ???
//		printf("\n\n Processo Iterativo no Tempo");
//		printf("\n Passo de tempo = %d \n",cont_t);
		//solucao inicial
		for(i = 0; i < neq; i++)
        {
            sol_ant[i] = sol[i];
            F_ant[i] = F[i];
        }

		montaSistemaGlobal_CN(nNosEl, A,AA,F,ID,IEN,COORD,BV,neq,nel,alpha,delta_t, ND, HN, tempo);//Galerkin + CN
		//Multiplica Matriz AA = (M - 0.5 dt K) por sol_ant e poe o resultado em P
		matvec(nNosEl, sol_ant,P,neq,nel,AA,ID,IEN);
		//calcula delta_t*(F[i] + F[i]) e acrescenta ao vetor P
		for (i = 0; i < neq; i++)
			P[i] = P[i] + alpha*delta_t*(F_ant[i] + F[i]);
		//Solução do sistema A*sol = P
		// solucaoGmres(A,P,sol,nel,neq,tol,kmax,lmax,ID,IEN);
		solucaoLCD(nNosEl, A,P,sol,nel,neq,tol,kmax,lmax,ID,IEN);
		tempo = tempo + delta_t;
		///////////////
		solG = criaVetorD(nnos);
		j = 0;
		for(i = 0; i < nnos; i++){
			if(ID[i] != 0){
				solG[i] = sol[j];
				j++;
			}
			else
				solG[i] = BV[i];
		}
/*
		sprintf(str, "galerkin%d.txt",cont_t);
		oout = fopen(str,"w");
		j=0;
		for(i = 0; i < nnos; i++){
			fprintf(oout,"%2.13lf  %2.13lf  %lf\n",COORD[i][0],COORD[i][1],solG[i]);
			if (j==nosx-1)//numero de nos
			{
				fprintf(oout,"\n");
				j=0;
			}
			else j++;

		}
		fclose(oout);
*/
		gravar(COORD, solG, nnos, nosx, nosy, cont_t);

		free(solG);

//		scriptFigGnuplot(nosx,nosy,limInfx,limSupx,limInfy,limSupy,cont_t);
//		scriptLatex(cont_t);
		//////////////
	}while(tempo < tfim);
	//FINAL DO PROCESSO ITERATIVO NO  TEMPO - ALGORITMO PREDITOR MULTICORRETOR



	tFinal = clock();
	tGasto = ((double)(tFinal - tInicio)/(CLOCKS_PER_SEC ));
	printf("\n\n ... TEMPO = %lf segundos\n",tGasto);

	printf("\n\n ... neq = %d\n",neq);
	printf("\n\n ... kmax = %d\n",kmax);

	printf("\n\n ... Quantidade de Iteracoes GMRES = %d\n",ITER_GMRES);
	printf("\n\n ... Quantidade de Iteracoes LCD = %d\n",ITER_LCD);



	solG = criaVetorD(nnos);
	j = 0;
	for(i = 0; i < nnos; i++){
		if(ID[i] != 0){
			solG[i] = sol[j];
			j++;
		}
		else
			solG[i] = BV[i];
	}

	free(sol);
	free(sol_ant);

	printf("\n\n ... Imprimindo resultados no arquivo saida");
/*
	oout = fopen("galerkin.txt","w");
	for(i = 0; i < nnos; i++)
		fprintf(oout,"%2.13lf  %2.13lf  %lf\n",COORD[i][0],COORD[i][1],solG[i]);
	fclose(oout);
	free(solG);

	scriptFigGnuplot(nosx,nosy,limInfx,limSupx,limInfy,limSupy);
*/

	//****************************************************************************
	// ...             liberando espaco de memoria utilizados
	//****************************************************************************
	printf("\n\n ... Liberando Espaco de Memoria\n");
	COORD = destruirMatrizD(COORD,nnos,2);
	IEN = destruirMatrizI(IEN,nel,nNosEl);
	free(BV);
	free(ID);
	free(F);
	free(P);
	//****************************************************************************

	//system("figSolucaoDV.plt");
	printf("\n\n");
//	system("pause");
//	system("gnuplot figsolucao.plt");
	return 0;
}
