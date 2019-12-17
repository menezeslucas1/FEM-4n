
#include<stdio.h>
#include<stdlib.h>
#include "matriz.h"
//===================================================================
//gera malha retangular
void gerarMalha(int nosx, int nosy, int tipoMalha, double limInfx, double limSupx, double limInfy, double limSupy, int **ien, double **coord){
	int i, j, k, e, no1, no2, no3, no4;
	int nel;
	double hx, hy;
	//FILE *f;

//	nel = ((nosx - 1)*(nosy - 1))*2;
	nel = ((nosx - 1)*(nosy - 1));

	hx = (limSupx - limInfx)/(nosx - 1.0);
	hy = (limSupy - limInfy)/(nosy - 1.0);

	// f = fopen("malha.txt","w");
	// fprintf(f,"%d  %d  %d  %d\n",nnos,nel,nosx,nosy);
	// fprintf(f,"\n");

	//malha inclinada para direita
	if(tipoMalha == 0){
		//numeracao global dos nohs
		no1 = 0;
		no2 = 0;
		no3 = 0;

		j = 0;
//		k = 2*nosx - 2; //numero de elemento em x
		k = nosx -1; //numero de elemento em x
		for(e = 1; e <= nel; e++){
			j++;
//			if(e%2 == 1){ //elemento eh impar
			no1++;
			no2 = no1 + 1;
			no3 = no2 + nosx;
			no4 = no1 + nosx;
//			}
/*
			else{ //elemento eh par
				no2 = no3;
				no3 = no2 - 1;
			}
*/
			//fprintf(f,"%d  %d  %d  %d\n",e,no1,no2,no3);
			ien[e - 1][0] = no1;
			ien[e - 1][1] = no2;
			ien[e - 1][2] = no3;
			ien[e - 1][3] = no4;

			if(j == k){
				no1++;
				j = 0;
			}
		}
		// fprintf(f,"\n");
/*		for (int e=1; e<=nel; e++){
			NLF[e-1][0]=0;
//			NLF[e-1][1]=0;
		}

		int i,j;
		for (i =1,j=1; i<=nosx-1;i++,j++){
			NLF[i*2-1-1][0]=j;
			NLF[i*2-1-1][1]=0;
		}
		NLF[--i*2-1][1]=++j;
		for(i=2; i<=nosy-1; i++,j++){
			NLF[i*(nosx-1)*2-1-1][0]=j;
			NLF[i*(nosx-1)*2-1-1][1]=0;
		}
		for (i=1; i<=nosy-1;i++,j++){
//			printf("nel=%d aux=%d\n",nel, (i-1)*(nosx-1)*2+2);
//			printf("i=%d\n",i);
			NLF[(i-1)*(nosx-1)*2+2-1][0]=j;
			NLF[(i-1)*(nosx-1)*2+2-1][1]=0;
		}
		NLF[(--i-1)*(nosx-1)*2+2-1][1]=++j;
		for(i=2; i<=nosx-1; i++,j++){
//			printf("nel=%d aux=%d\n",nel, (nosy-2)*(nosx-1)*2 + i*2);
//			printf("i=%d\n",i);
			NLF[(nosy-2)*(nosx-1)*2 + i*2-1][0]=j;
			NLF[(nosy-2)*(nosx-1)*2 + i*2-1][1]=0;
		}
*/
	}
	else{ //Malha do tipo 1 *** inclinacao para esquerda
		no1 = 1;
		no2 = 0;
		no3 = 0;

		j = 0;
		k = 2*nosx - 2; //numero de elemento em x
		for(e = 1; e <= nel; e++){
			j++;
			if(e%2 == 1){ //elemento eh impar
				no1 = no1;
				no2 = no1 + 1;
				no3 = no2 + nosx - 1;
			}
			else{ //elemento eh par
				no1 = no2;
				no2 = no3 + 1;
				no3 = no2 - 1;
			}
			//fprintf(f,"%d  %d  %d  %d\n",e,no1,no2,no3);

			ien[e - 1][0] = no1;
			ien[e - 1][1] = no2;
			ien[e - 1][2] = no3;

			if(j == k){
				no1++;
				j = 0;
			}
		}
	}

	//coordenadas dos nos globais
	k = 1;
	for(i = 1; i <= nosy; i++){
		for(j = 1; j <= nosx; j++){
			//fprintf(f,"%d  %f  %f\n",k,limInfx + (j-1)*hx,limInfy + (i-1)*hy);
			coord[k-1][0] = limInfx + (j-1)*hx;
			coord[k-1][1] = limInfy + (i-1)*hy;
			k++;
		}
	}
	// fclose(f);
/*
	j=1;
	for (i =1; i<=nosx-1; i++,j++){
		front[j].v[0]=i;
		front[j].v[1]=i+1;
	}
	for (i=1; i<=nosy-1; i++,j++){
		front[j].v[0]=i*nosx;
		front[j].v[1]=(i+1)*nosx;
	}
	for(i=1; i<=nosy-1; i++,j++){
		front[j].v[0]=(i-1)*nosx+1;
		front[j].v[1]=i*nosx+1;
	}
	for(i=1; i<=nosx; i++, j++){
		front[j].v[0]=(nosy-1)*nosx+i;
		front[j].v[0]=(nosy-1)*nosx+i+1;
	}
*/
}

//==============================================================================
//preenche o vetor  ID
int condFront(int *id, int *nd,int **hn, int nnos, int nosx, int nosy, int opcao){
	int k, i, j, neq;
	//CONDICOES DE FRONTEIRAS

	//noh = 0 condicao de fronteira prescrita
	//noh = 1 fronteira nao prescrita

	neq = 0;
	printf("condfront opcao=%d\n",opcao);
	switch(opcao){
	case 1: //fronteira prescrita em todo o contorno
		//printf("\n ... OPCAO 3 - fronteira prescrita em todo o contorno");
		k = 1;
		for(j = 1; j <= nosy; j++){
			for(i = 1; i <= nosx; i++){
				if((j==1)||(i==1)||(i==nosx)||(j==nosy)){
					id[k-1] = 0;
					nd[k-1] = 0;
//					front[k-1]=p;
//					p++;
				}
				else{
					neq++;

					id[k-1] = neq;
				}
				k++;
			}
		}
		break;
	case 2://fronteira com um lado de Neumann ( baixo )
		k = 1;
		for(j = 1; j <= nosy; j++){
			for(i = 1; i <= nosx; i++){
				if((i==1)||(i==nosx)||(j==nosy)){
					id[k-1] = 0;
					nd[k-1] = 0;
//					front[k-1]=p;
//					p++;
				}
				else{
					neq++;
					id[k-1] = neq;
					if (j==1){
						nd[k-1] = 1;//nó de Neumann
//						front[k-1]=p;
//						p++;
					}
					else
						nd[k-1] = 0;
				}
				k++;
			}
		}
		hn[0][0]=-1;
		hn[0][1]=-1;
		hn[1][0]= 2;
		hn[1][1]=-1;
		for (k=2; k<=nosx-3; k++){
            hn[k][0]=k-1;
            hn[k][1]=k+1;
        }
        hn[nosx-2][0]=nosx-3;
        hn[nosx-2][1]=-1;

        hn[nosx-1][0]=-1;
        hn[nosx-1][1]=-1;
		for (k=nosx; k<nnos; k++){
            hn[k][0]=-1;
            hn[k][1]=-1;
        }
		break;
	}
	return neq;
}


//===================================================================
//a contagem é feita no sentido anti-horario
// retorna em qual fronteira esta o elemento
// retorna -1 se elemento nao estah na fronteira
// retorna 1 se elemento esta na fronteira 1
// retorna 2 se elemento esta na fronteira 2
// retorna 3 se elemento esta na fronteira 3
// retorna 4 se elemento esta na fronteira 4
// retorna 5 se elemento esta no corner 5
// retorna 6 se elemento esta no corner 6
//limInfx,limSupx,limInfy,limSupy
int fronteiraSaida(int e, int **ien, double **coordenada){
	int front, noLocal, noGlobal;
	double *x, *y, F1, F2, F3, F4;
	x = criaVetorD(3); //x = [x0, x1, x2]
	y = criaVetorD(3); //y = [y0, y1, y2]

	for(noLocal = 0; noLocal < 3; noLocal++){
		noGlobal = ien[e-1][noLocal];
		x[noLocal] = coordenada[noGlobal - 1][0]; //coordenada x do noh global
		y[noLocal] = coordenada[noGlobal - 1][1]; //coordenada y do noh global

	}

	F1 = 0.0;
	F2 = 0.0;
	F3 = 1.0;
	F4 = 1.0;

	if((y[0] == F1)&&(y[1] == F1)){ //fronteira1
		if(x[1] != F4)     //eixo x = fronteira 1
			front = 1;
		else
			front = 5;      //corner 5
	}
	else{//1o else
		if((x[0] == F2)&&(x[2] == F2)){ //fronteira2
			if(y[2] != F3)     //eixo y = fronteira 2
				front = 2;
			else
				front = 6;      //corner 6
		}
		else{//2o else
			if((y[1] == F3)&&(y[2] == F3)) //fronteira3
				front = 3;
			else{//3o else
				if((x[1] == F4)&&(x[2] == F4)) //fronteira4
					front = 4;
				else //elemento nao estah na fronteira
					front = -1;
			}//fechando 3o else
		}//fechando 2o else
	}//fechando 1o else

	if(front == 0)
		printf("\nElemento = %d - Fronteira = %d",e,front);

	free(x); free(y);
	return front;

}
