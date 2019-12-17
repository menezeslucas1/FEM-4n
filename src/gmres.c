
//Galerkin Method
//===================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matriz.h"
#include <gsl/gsl_cblas.h>

extern int ITER_GMRES;

//=======================================================================================================

void dclear(double *a , int n){
	int i;

	for(i = 0; i < n; i++)
		a[i] = 0.0;
}

//========================================================================================================

void dclear_matrix(double **a , int m, int n){
	int i, j;

	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			a[i][j] = 0.0;

}

//========================================================================================================

void dclear_matrix_int(int **a , int m, int n){
	int i, j;

	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			a[i][j] = 0;
}
//=========================================================================================================

//===================================================================
void matvec(int nNosEl, double *p, double *ap, int neq, int nel, double ***s, int *id, int **ien){
	int i, j, e, neqv[nNosEl];
	double pv[nNosEl], apv[nNosEl];

	dclear(ap,neq);

	//loop over element blocks
	for(e = 1; e <= nel; e++){
		//gather global components of p
		for(i = 0; i < nNosEl; i++){
			neqv[i] = id[ien[e - 1][i] - 1]; //equacao associada ao noh local i  + 1 = 1,2,3
			//if(e == 11)
			// printf("\n i = %d e eq = %d",i,neqv[i]);
		}
		for(i = 0; i < nNosEl; i++){
			if(neqv[i] == 0)
				pv[i] = 0.0;
			else
				pv[i] = p[neqv[i]-1];
		}
		//ebe matrix-vector multiplication for all elements in this block
		for(i = 0; i < nNosEl; i++){
			apv[i] = 0.0;
			for(j = 0; j < nNosEl; j++)
				apv[i] += s[e-1][i][j]*pv[j];
		}
		//scatter and add to global 'ap'
		for(i = 0; i < nNosEl; i++){
			if(neqv[i] != 0)
				ap[neqv[i]-1] += apv[i];
		}
	}
}

//======================================================================================================
void solucaoGmres(int nNosEl, double ***A, double *b, double *x, int nel, int neq, double etol, int kmax, int lmax, int *id, int **ien){
	int i, j, k, l, cont;
	double xnormb, eps, hbji, hbjk, r, rho, soma;
	double **u, **Hbar, *ce, *es, *yip, *e;

	u = criaMatrizD(kmax,neq);
	Hbar = criaMatrizD(kmax,kmax);
	e = criaVetorD(kmax);
	ce = criaVetorD(kmax);
	es = criaVetorD(kmax);
	yip = criaVetorD(kmax);

	//compute ||b||_2
	xnormb = cblas_ddot(neq,b,1,b,1);
	//imprimirVetorD(b,neq);
	xnormb = sqrt(xnormb);

	//clear u
	dclear_matrix(u,kmax,neq);

	//clear x  (x <- 0)
	dclear(x,neq);

	//computa eps = eps_tol*b
	eps = etol*xnormb;

	l = 0;
	cont = 0;

	//begin GMRES cycles
	do{
		//printf("\n\n =========================== ITERACAO l =  %d ==============================",l);
		i = 0;
		//ebe matrix-vector multiplication ... ui = Ax
		matvec(nNosEl,x,u[i],neq,nel,A,id,ien);
		cblas_daxpy(neq,-1.0,b,1,u[i],1);   //u_i <--  u_i - b
		cblas_dscal(neq,-1,u[i],1);         //u_i <-- -u_i = (b - u_i)

		//calcula e_i = ||u_i||
		e[i] = cblas_ddot(neq,u[i],1,u[i],1);
		e[i] = sqrt(e[i]);

		//u_i = u_i/e_i
		cblas_dscal(neq,1.0/e[i],u[i],1);

		rho = e[i];

		while((rho > eps)&&(i < kmax - 1)){
			// printf("\n\n Iteracao i =  %d",i);
			cont = cont + 1;
			//ebe matrix-vector multiplication ... uj = u_i+1 = Aui
			matvec(nNosEl, u[i],u[i+1],neq,nel,A,id,ien);
			//modified gram-schimidt reorthogonalization
			for(j = 0; j <= i; j++){
				//produto interno Hji = (u_i+1, u_j)
				Hbar[j][i] = cblas_ddot(neq,u[j],1,u[i+1],1);  //ui = u_i+1
				cblas_daxpy(neq,-Hbar[j][i],u[j],1,u[i+1],1);   //u1 = u_i+1 <-- u_i+1 - hbar*u0
			}//end_for j
			//produto interno H_i+1_i = (u_i+1, u_i+1)
			Hbar[i+1][i] = cblas_ddot(neq,u[i+1],1,u[i+1],1);
			Hbar[i+1][i] = sqrt(Hbar[i+1][i]);

			//calcula u_i+1 = u_i+1/H_i+1,i
			cblas_dscal(neq,1.0/Hbar[i+1][i],u[i+1],1);   //u_i+1 = ui  <-- ui/h_i+1,i
			//Q-R algorithm
			//if(iter > 1){
			for(j = 0; j <= i-1; j++){
				hbji =  ce[j]*Hbar[j][i] + es[j]*Hbar[j+1][i];
				hbjk = -es[j]*Hbar[j][i] + ce[j]*Hbar[j+1][i];
				Hbar[j][i]   = hbji;
				Hbar[j+1][i] = hbjk;
			}//j

			r = Hbar[i][i]*Hbar[i][i] +  Hbar[i+1][i]*Hbar[i+1][i];
			r = sqrt(r);
			ce[i] = Hbar[i][i]/r;
			es[i] = Hbar[i+1][i]/r;
			Hbar[i][i] = r;
			Hbar[i+1][i] = 0.0;
			e[i+1] = -es[i]*e[i];
			e[i]   =  ce[i]*e[i];
			rho = fabs(e[i+1]);
			//printf("\n residuo  = %.10lf",rho);
			i = i + 1;
		}//end_while rho i
		i = i - 1;
		yip[i] = e[i]/Hbar[i][i];
		for(j = i - 1; j >= 0; j--){
			soma = 0.0;
			for(k = j + 1; k <= i; k++)
				soma = soma + Hbar[j][k]*yip[k];
			yip[j] = (e[j] - soma)/Hbar[j][j];
		}
		for(j = 0; j <= i; j++)
			for(k = 0; k < neq; k++)
				x[k] = x[k] + u[j][k]*yip[j];

		l = l + 1;

	}while((rho > eps)&&(l<lmax));

	//printf("\n\n residuo  = %.10lf *** eps = %f *** tol = %f *** xnormb = %f *** l = %d *** lmax = %d",rho,eps,etol,xnormb,l,lmax);
	ITER_GMRES =  ITER_GMRES + cont;
	u = destruirMatrizD(u,kmax,neq);
	Hbar = destruirMatrizD(Hbar,kmax,kmax);

	free(e);
	free(ce);
	free(es);
	free(yip);
}
