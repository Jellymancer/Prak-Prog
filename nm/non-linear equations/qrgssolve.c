#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<stdio.h>
#include<gsl/gsl_blas.h>


void qrgssolve(gsl_matrix* Qqt,  gsl_matrix* Rr,gsl_vector* b,gsl_vector* x){

// Solving the system R*x=Qt*b by back substitution


// Define c=Qt*b
gsl_vector* c=gsl_vector_calloc(x->size);

gsl_blas_dgemv(CblasNoTrans,1.0,Qqt,b,0.0,c);


for(int i=c->size-1; i>=0; i--){
	double s=gsl_vector_get(c,i);
	for(int k=i+1;k< c->size; k++)
		s-=gsl_matrix_get(Rr,i,k)*gsl_vector_get(c,k);
		gsl_vector_set(c,i,s/gsl_matrix_get(Rr,i,i));}

gsl_vector_memcpy(x,c);
gsl_vector_free(c);
}
