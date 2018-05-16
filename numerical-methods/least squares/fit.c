#include "fit.h"
#include "matrixprep.h"
#include<gsl/gsl_linalg.h>

void lsfit(gsl_matrix* M, gsl_vector* b, gsl_vector* c, gsl_matrix* Sigma){
int n=M->size1; int m=M->size2;

gsl_matrix* Q=gsl_matrix_calloc(n,m); // allocating memory
gsl_matrix* R=gsl_matrix_calloc(m,m);
gsl_matrix* Sigma2=gsl_matrix_calloc(m,m);
gsl_matrix* SigmaQ=gsl_matrix_calloc(m,m);


gsl_matrix_memcpy(Q,M);
// solving c=R⁻¹Q^Tb using QR decomposition.
qrgsdecomp(Q,R);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,R,R,0,Sigma2);
qrgssolve(Q,R,b,c);


gsl_vector* q=gsl_vector_calloc(m);
gsl_vector* x=gsl_vector_calloc(m);
gsl_matrix_memcpy(SigmaQ,Sigma2);
qrgsdecomp(SigmaQ,R);

for(int i=0;i< m;i++){ //Calculating the covariance matrix
        gsl_vector_set(q,i,1);
        qrgssolve(SigmaQ,R,q,x);
        gsl_matrix_set_col(Sigma,i,x);
        gsl_vector_set_zero(q);
        gsl_vector_set_zero(x);
}


gsl_matrix_free(Q);
gsl_matrix_free(R);
gsl_matrix_free(Sigma2);
gsl_matrix_free(SigmaQ);
gsl_vector_free(q);
gsl_vector_free(x);
}

//Rest of the function are used for QR decomp and solving the decomposed equation.
void qrgsdecomp(gsl_matrix *E,gsl_matrix *W){

int s=E->size2; //Number of columns in matrix E
for(int i=0;i<s;i++){
	gsl_vector_view col=gsl_matrix_column(E,i);
	double Rii = gsl_blas_dnrm2(&col.vector);
	gsl_matrix_set(W,i,i,Rii);
	gsl_vector_scale(&col.vector,1/Rii);

	for(int j=i+1;j<s;j++){
		gsl_vector_view col2= gsl_matrix_column(E,j);
		double Rij = 0;
		gsl_blas_ddot(&col.vector,&col2.vector,&Rij);
		gsl_blas_daxpy(-Rij,&col.vector,&col2.vector);
		gsl_matrix_set(W,i,j,Rij);
		}
	}

}


void qrgssolve(gsl_matrix* Q,  gsl_matrix* R,gsl_vector* b,gsl_vector* c){
gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,c);
for(int i=c->size-1; i>=0; i--){
	double s=gsl_vector_get(c,i);
	for(int k=i+1;k< c->size; k++)
		s-=gsl_matrix_get(R,i,k)*gsl_vector_get(c,k);
		gsl_vector_set(c,i,s/gsl_matrix_get(R,i,i));}

}



