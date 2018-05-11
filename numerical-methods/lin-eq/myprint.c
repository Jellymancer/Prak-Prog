#include<stdlib.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<stdio.h>


void printm(gsl_matrix *T){
for(int c=0;c< T->size1;c++){
	for(int r=0;r< T->size2;r++){
		printf("%7.3f  ",gsl_matrix_get(T,c,r));
}
printf("\n");
}
printf("\n");
}

void printv(gsl_vector *v){
	for(int r=0;r< v->size;r++){
		printf("%7.3f",gsl_vector_get(v,r));
		printf("\n");}
}
