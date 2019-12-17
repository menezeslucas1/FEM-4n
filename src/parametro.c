
#include<stdio.h>
#include<math.h>

void solucaoInicial(double *sol, int *ID, double **COORD, int nnos){
	int i, eq;

	for(i = 0; i < nnos; i++){
		eq = ID[i];   //noh global associado a noLocal
		if(eq!=0){
			sol[eq - 1] = 0.0;
		}
	}
}

//===================================================================
double funcaoFonte(double x, double y){
	double valor = exp(-100.0*((x-0.3)*(x-0.3) + (y-0.3)*(y-0.3)));
	return valor; //fonte variável
//	return 1.0; //fonte constante
}

//===================================================================
double difusao(){
//	return pow(10.0,-6.0);
	return pow(10.0,-1.0);
}

//===================================================================
double reacao(){
	return 0.0;
}

//===================================================================
double velocx(double x, double y, double t){
//	return 0; //constante
//	return 1;
	return -1*exp(10.0*t)*(y-0.5);// variável
}

//===================================================================
double velocy(double x, double y, double t){
//	return -1;
//	return 0.5;
//	return 0;
	return 1*exp(10.0*t)*(x-0.5);
}

//====================================================================
//Parametro Cb para calculo da solucao inicial
double parametroCb(){
	return 0.4;
}

//====================================================================
//kmax gmres
int kmax_gmres(){
	return 50;
}

//====================================================================
//lmax gmres
int lmax_gmres(){
	return 50;
}

//====================================================================
//tolerancia gmres
double tol_gmres(){
	return pow(10.0,-3.0);;
}

//====================================================================
//kkmax metodo DV
int kkmax_metodoDV(){
	return 50;
}


//====================================================================
//preenche o vetor BV(condição de fronteira)
void lerCondFront(double *bv, int *id, int nosx, int nosy){
	int i;
	int nnos = nosx*nosy;

	for(i = 1; i <= nnos; i++){
		bv[i-1]= 0.0;    //valor de fronteira de cada noh eh zero
		if(id[i-1] == 0){ //noh i prescrito

			//condicao contorno - campo de velocidade (1,1) ou (0.5, 1)

			//Paredao - Exemplo 2
			/*if((i >= nosx/3 + 1)&&(i <= nosx))
				bv[i-1] = 1.0;
			else
				bv[i-1] = 0.0;
			*/

			//Problema das 02 rampas
			//if((i >= nosx/2 + 1)&&(i <= nosx))
			//   bv[i-1] = 1.0;
			//else
			//   bv[i-1] = 0.0;
			//Problema de 01 rampa - condicao contorno = 0 - Exemplo 3
			bv[i-1] = 0.0;
		}
	}
}

double gn(double x,double y){
	return 0;
//	return +0.1;
//	return 1;
//	return -1;
}
