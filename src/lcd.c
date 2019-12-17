#include<stdio.h>
#include<math.h>
#include "gmres.h"
#include "matriz.h"
#include <gsl/gsl_blas.h>

extern int ITER_LCD;

//======================================================================================================

//copia vetor u em v
void copy_vector(double *u, double *v, int n){
	int i;

	for(i = 0; i < n; i++)
		v[i] = u[i];
}

//======================================================================================================

void solucaoLCD(int nNosEl, double ***A, double *b, double *x, int nel, int neq, double etol, int kmax, int lmax, int *id, int **ien){
	int i, j, l, cont;
	double eps, rho, alpha, beta;
	double **P, **Q;
	double *res;
	P = criaMatrizD(kmax,neq);
	Q = criaMatrizD(kmax,neq);
//printf("oi\n");
	res = criaVetorD(neq);
	dclear(x,neq); // x <-- 0
	copy_vector(b,res,neq); //res <--- b
	rho = cblas_ddot(neq,res,1,res,1); //rho <-- ||res||_2
	rho = sqrt(rho);
	eps = etol*rho;   //eps = eps_tol*||res||_2
	l = 0;
	cont = 0;
	copy_vector(res,P[0],neq); //p[0] <--- res

	//begin LCD cycles
	do{
		//printf("\n\n =========================== ITERACAO l =  %d ==============================",l);
		//ebe matrix-vector multiplication ... ql = A pl
		matvec(nNosEl, P[0],Q[0],neq,nel,A,id,ien);
		i = 0;
		while((rho > eps)&&(i < kmax - 1)){
			cont = cont + 1;
			alpha = cblas_ddot(neq,P[i],1,res,1)/cblas_ddot(neq,P[i],1,Q[i],1);
			cblas_daxpy(neq,alpha,P[i],1,x,1);          //x <--  x + alpha*pi
			cblas_daxpy(neq,-1.0*alpha,Q[i],1,res,1);   //res <--  res - alpha*qi
			copy_vector(res,P[i+1],neq);       //p[i+1] <-- res
			matvec(nNosEl, P[i+1],Q[i+1],neq,nel,A,id,ien);  //q[i+1] = A p[i+1]

			for(j = 0; j <= i; j++){
				beta = -1.0*(cblas_ddot(neq,P[j],1,Q[i+1],1)/cblas_ddot(neq,P[j],1,Q[j],1));
				cblas_daxpy(neq,beta,P[j],1,P[i+1],1);   //p[i+1] <-- p[i+1] + beta*p[j]
				cblas_daxpy(neq,beta,Q[j],1,Q[i+1],1);   //q[i+1] <-- q[i+1] + beta*q[j]
			}//end_for j

			//compute ||res||_2
			rho = cblas_ddot(neq,res,1,res,1);
			rho = sqrt(rho);
			//printf("\n residuo LCD  = %.10lf",rho);
			i = i + 1;
		}//end_while rho i
		//i = i - 1;
		copy_vector(P[i],P[0],neq);  //p[0] <-- p[kmax]
		l = l + 1;
	}while((rho > eps)&&(l<lmax));
	//printf("\n\n residuo  = %.10lf *** eps = %f *** tol = %f *** xnormb = %f *** l = %d *** lmax = %d",rho,eps,etol,xnormb,l,lmax);
	ITER_LCD =  ITER_LCD + cont;
	P = destruirMatrizD(P,kmax,neq);
	Q = destruirMatrizD(Q,kmax,neq);
	free(res);
}
