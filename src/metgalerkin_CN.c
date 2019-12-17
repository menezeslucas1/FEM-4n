
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "matriz.h"
#include "parametro.h"
#include <math.h>

//==============================================================================
void matrizVetorLocal_CN(int nNosEl, int e, double alpha, double delta_t, double **E, double **EE, double *Fe, int *id, int **ien, double **coordenada, double *bv, int *ND, int **HN, double t){
	int i, j, noLocal, noGlobal;
	double *x, *y, **M, **K;
	double difus, reac, vx, vy, h, delta, modv, csi, peclet;
	double a, b;

	x = criaVetorD(nNosEl); //x = [x0, x1, x2]
	y = criaVetorD(nNosEl); //y = [y0, y1, y2]

	//Preenchendo vetores x e y - guardam as coordenadas dos nohs globais do elemento e
	for(noLocal = 0; noLocal < nNosEl; noLocal++){
		noGlobal = ien[e-1][noLocal];   //noh global associado a noLocal
		x[noLocal] = coordenada[noGlobal - 1][0]; //coordenada x do noh global
		y[noLocal] = coordenada[noGlobal - 1][1]; //coordenada y do noh global
	}
	a = (x[1] - x[0])/2;
	b = (y[2] - y[0])/2;
	h = a*2;
//	printf("a=%f, b=%f\n",a,b);

	difus = difusao();
	reac = reacao();
	vx = velocx((x[0]+x[1]+x[2]+x[3])/4.0, (y[0]+y[1]+y[2]+y[3])/4.0, t);
	vy = velocy((x[0]+x[1]+x[2]+x[3])/4.0, (y[0]+y[1]+y[2]+y[3])/4.0, t);
//	printf("vx=%f vy=%f\n", vx,vy);
	modv = sqrt(vx*vx + vy*vy);
	peclet = modv*h/(2*difus);
	csi = (1-1/peclet);
	csi = (0 >= csi)?0:csi;
    delta = h*csi/(2*modv);

	M = criaMatrizD(nNosEl,nNosEl);
	K = criaMatrizD(nNosEl,nNosEl);

	//Matriz M
	M[0][0] = -b * delta * vx / 0.3e1 + 0.4e1 / 0.9e1 * a * b - a * delta * vy / 0.3e1;
	M[0][1] = -b * delta * vx / 0.3e1 + 0.2e1 / 0.9e1 * a * b - a * delta * vy / 0.6e1;
	M[0][2] = -b * delta * vx / 0.6e1 + a * b / 0.9e1         - a * delta * vy / 0.6e1;
	M[0][3] = -b * delta * vx / 0.6e1 + 0.2e1 / 0.9e1 * a * b - a * delta * vy / 0.3e1;
	M[1][0] =  b * delta * vx / 0.3e1 + 0.2e1 / 0.9e1 * a * b - a * delta * vy / 0.6e1;
	M[1][1] =  b * delta * vx / 0.3e1 + 0.4e1 / 0.9e1 * a * b - a * delta * vy / 0.3e1;
	M[1][2] =  b * delta * vx / 0.6e1 + 0.2e1 / 0.9e1 * a * b - a * delta * vy / 0.3e1;
	M[1][3] =  b * delta * vx / 0.6e1 + a * b / 0.9e1         - a * delta * vy / 0.6e1;
	M[2][0] =  b * delta * vx / 0.6e1 + a * b / 0.9e1         + a * delta * vy / 0.6e1;
	M[2][1] =  b * delta * vx / 0.6e1 + 0.2e1 / 0.9e1 * a * b + a * delta * vy / 0.3e1;
	M[2][2] =  b * delta * vx / 0.3e1 + 0.4e1 / 0.9e1 * a * b + a * delta * vy / 0.3e1;
	M[2][3] =  b * delta * vx / 0.3e1 + 0.2e1 / 0.9e1 * a * b + a * delta * vy / 0.6e1;
	M[3][0] = -b * delta * vx / 0.6e1 + 0.2e1 / 0.9e1 * a * b + a * delta * vy / 0.3e1;
	M[3][1] = -b * delta * vx / 0.6e1 + a * b / 0.9e1         + a * delta * vy / 0.6e1;
	M[3][2] = -b * delta * vx / 0.3e1 + 0.2e1 / 0.9e1 * a * b + a * delta * vy / 0.6e1;
	M[3][3] = -b * delta * vx / 0.3e1 + 0.4e1 / 0.9e1 * a * b + a * delta * vy / 0.3e1;
	//******************************************************************************************************
	//Matriz K
	double a2 = a*a;
	double b2 = b*b;
	double ab = a*b;
	double vx2 = vx * vx;
	double vy2 = vy * vy;
	double dr = delta * reac;
	double a2bdrvy = a2 * b * dr * vy;
	double ab2drvx = a * b2 * dr * vx;
	double a2b2r = a2 * b2 * reac;
	double a2dvy2 = a2 * delta * vy2;
	double b2dvx2 = b2 * delta * vx2;
	double a2bvy = a2 * b * vy;
	double ab2vx = a * b2 * vx;
	double a2d = a2 * difus;
	double b2d = b2 * difus;
	double abdvxvy = ab * delta * vx * vy;
	K[0][0] = (double) (-6 * a2bdrvy - 6 * ab2drvx + 8 * a2b2r + 6 * a2dvy2 + 6 * b2dvx2 - 6 * a2bvy - 6 * ab2vx + 6 * a2d + 6 * b2d + 9 * abdvxvy) / (0.18e2 * ab);
	K[0][1] = (double) (-3 * a2bdrvy - 6 * ab2drvx + 4 * a2b2r + 3 * a2dvy2 - 6 * b2dvx2 - 3 * a2bvy + 6 * ab2vx + 3 * a2d - 6 * b2d              ) / (0.18e2 * ab);
	K[0][2] = (double) (-3 * a2bdrvy - 3 * ab2drvx + 2 * a2b2r - 3 * a2dvy2 - 3 * b2dvx2 + 3 * a2bvy + 3 * ab2vx - 3 * a2d - 3 * b2d - 9 * abdvxvy) / (0.18e2 * ab);
	K[0][3] = (double) (-6 * a2bdrvy - 3 * ab2drvx + 4 * a2b2r - 6 * a2dvy2 + 3 * b2dvx2 + 6 * a2bvy - 3 * ab2vx - 6 * a2d + 3 * b2d              ) / (0.18e2 * ab);
	K[1][0] = (double) (-3 * a2bdrvy + 6 * ab2drvx + 4 * a2b2r + 3 * a2dvy2 - 6 * b2dvx2 - 3 * a2bvy - 6 * ab2vx + 3 * a2d - 6 * b2d              ) / (0.18e2 * ab);
	K[1][1] = (double) (-6 * a2bdrvy + 6 * ab2drvx + 8 * a2b2r + 6 * a2dvy2 + 6 * b2dvx2 - 6 * a2bvy + 6 * ab2vx + 6 * a2d + 6 * b2d - 9 * abdvxvy) / (0.18e2 * ab);
	K[1][2] = (double) (-6 * a2bdrvy + 3 * ab2drvx + 4 * a2b2r - 6 * a2dvy2 + 3 * b2dvx2 + 6 * a2bvy + 3 * ab2vx - 6 * a2d + 3 * b2d              ) / (0.18e2 * ab);
	K[1][3] = (double) (-3 * a2bdrvy + 3 * ab2drvx + 2 * a2b2r - 3 * a2dvy2 - 3 * b2dvx2 + 3 * a2bvy - 3 * ab2vx - 3 * a2d - 3 * b2d + 9 * abdvxvy) / (0.18e2 * ab);
	K[2][0] = (double) ( 3 * a2bdrvy + 3 * ab2drvx + 2 * a2b2r - 3 * a2dvy2 - 3 * b2dvx2 - 3 * a2bvy - 3 * ab2vx - 3 * a2d - 3 * b2d - 9 * abdvxvy) / (0.18e2 * ab);
	K[2][1] = (double) ( 6 * a2bdrvy + 3 * ab2drvx + 4 * a2b2r - 6 * a2dvy2 + 3 * b2dvx2 - 6 * a2bvy + 3 * ab2vx - 6 * a2d + 3 * b2d              ) / (0.18e2 * ab);
	K[2][2] = (double) ( 6 * a2bdrvy + 6 * ab2drvx + 8 * a2b2r + 6 * a2dvy2 + 6 * b2dvx2 + 6 * a2bvy + 6 * ab2vx + 6 * a2d + 6 * b2d + 9 * abdvxvy) / (0.18e2 * ab);
	K[2][3] = (double) ( 3 * a2bdrvy + 6 * ab2drvx + 4 * a2b2r + 3 * a2dvy2 - 6 * b2dvx2 + 3 * a2bvy - 6 * ab2vx + 3 * a2d - 6 * b2d              ) / (0.18e2 * ab);
	K[3][0] = (double) ( 6 * a2bdrvy - 3 * ab2drvx + 4 * a2b2r - 6 * a2dvy2 + 3 * b2dvx2 - 6 * a2bvy - 3 * ab2vx - 6 * a2d + 3 * b2d              ) / (0.18e2 * ab);
	K[3][1] = (double) ( 3 * a2bdrvy - 3 * ab2drvx + 2 * a2b2r - 3 * a2dvy2 - 3 * b2dvx2 - 3 * a2bvy + 3 * ab2vx - 3 * a2d - 3 * b2d + 9 * abdvxvy) / (0.18e2 * ab);
	K[3][2] = (double) ( 3 * a2bdrvy - 6 * ab2drvx + 4 * a2b2r + 3 * a2dvy2 - 6 * b2dvx2 + 3 * a2bvy + 6 * ab2vx + 3 * a2d - 6 * b2d              ) / (0.18e2 * ab);
	K[3][3] = (double) ( 6 * a2bdrvy - 6 * ab2drvx + 8 * a2b2r + 6 * a2dvy2 + 6 * b2dvx2 + 6 * a2bvy - 6 * ab2vx + 6 * a2d + 6 * b2d - 9 * abdvxvy) / (0.18e2 * ab);

	//Vetor Fe
	double f1 = funcaoFonte(coordenada[ien[e - 1][0] - 1][0],coordenada[ien[e - 1][0] - 1][1]);
	double f2 = funcaoFonte(coordenada[ien[e - 1][1] - 1][0],coordenada[ien[e - 1][1] - 1][1]);
	double f3 = funcaoFonte(coordenada[ien[e - 1][2] - 1][0],coordenada[ien[e - 1][2] - 1][1]);
	double f4 = funcaoFonte(coordenada[ien[e - 1][3] - 1][0],coordenada[ien[e - 1][3] - 1][1]);
	Fe[0] = (double) (a * b * (4 * f1 + 2 * f2 + f3 + 2 * f4)) / 0.9e1;
	Fe[1] = (double) (a * b * (2 * f1 + 4 * f2 + 2 * f3 + f4)) / 0.9e1;
	Fe[2] = (double) (a * b * (f1 + 2 * f2 + 4 * f3 + 2 * f4)) / 0.9e1;
	Fe[3] = (double) (a * b * (2 * f1 + f2 + 2 * f3 + 4 * f4)) / 0.9e1;

	//Matriz do elemento
	for(i = 0; i < nNosEl; i++)
		for(j = 0; j < nNosEl; j++){
			E[i][j] =  M[i][j] + (delta_t*K[i][j])/2.0;
			EE[i][j] =  M[i][j] - (delta_t*K[i][j])/2.0;
		}

	//condicoes de fronteira no vetor Fe
/*
	for(i = 0; i < nNosEl; i++){
		noGlobal = ien[e-1][i];  //noGlobal associado ao noLocal i
		if(id[noGlobal - 1] != 0) //se equacao do noGlobal eh diferente de zero
		{
//			Fe[i] = Fe[i] - (bv[ien[e-1][0] - 1]*K[i][0] + bv[ien[e-1][1] - 1]*K[i][1] + bv[ien[e-1][2] - 1]*Ae[i][2]);
			if (ND[noGlobal-1]==1)// nó na fronteira de Neumann
				for (j=0; j<nNosEl; j++)//procura aresta de Neumann no elemento
					if (ND[ien[e-1][j]]==1)//verifica se tem outro vértice de Neumann no elemento
						if(HN[noGlobal-1][0]==ien[e-1][j] || HN[noGlobal-1][1]==ien[e-1][j]){//verifica se esses vértices formam uma aresta de Neumann
							double x1,y1,x2,y2;
							x1 = coordenada[noGlobal-1][0];
							x2 = coordenada[ien[e-1][j]][0];
							y1 = coordenada[noGlobal-1][1];
							y2 = coordenada[ien[e-1][j]][1];
							double sx = x1 - x2;
							sx = sx*sx;
							double sy = y1 - y2;
							sy = sy*sy;
							double s = sqrt(sx + sy);
							double fgn = gn(coordenada[noGlobal-1][0], coordenada[noGlobal-1][1]) * s /2;
							Fe[i] += fgn;
						}
		}
	}
*/
	M = destruirMatrizD(M,nNosEl,nNosEl);
	K = destruirMatrizD(K,nNosEl,nNosEl);
	free(x); free(y);
}

//==============================================================================
void montaSistemaGlobal_CN(int nNosEl, double ***A, double ***AA, double *F, int *id, int **ien, double **coordenada, double *bv, int neq, int nel, double alpha, double delta_t, int *ND, int **HN, double t){
	double **E, **EE, *Fe;
	//  int i, e, eq, a, b;
	int e, eq, a, b;

	E = criaMatrizD(nNosEl,nNosEl); //matriz local
	EE = criaMatrizD(nNosEl,nNosEl); //matriz local
	Fe = criaVetorD(nNosEl);    //vetor forca local

	for (a =0; a<neq; a++)
	{
		F[a] = 0;
	}
	//visitando cada elemento da malha
	for(e = 1; e <= nel; e++){
		matrizVetorLocal_CN(nNosEl, e,alpha,delta_t,E,EE,Fe,id,ien,coordenada,bv, ND, HN, t);
		//Preenchendo vetor global F
		for(a = 0; a < nNosEl; a++){
			eq = id[ien[e-1][a] - 1];//equacao associada ao noh local a
			if(eq != 0){
				F[eq - 1] += Fe[a];
			}
		}
		//Preenchendo matrizes globais tridimensional A e AA
		for(a = 0; a < nNosEl; a++)
			for(b = 0; b < nNosEl; b++){
				A[e - 1][a][b] = E[a][b];
				AA[e - 1][a][b] = EE[a][b];
			}
	}

	E = destruirMatrizD(E,nNosEl,nNosEl);
	EE = destruirMatrizD(EE,nNosEl,nNosEl);
	free(Fe);
}
