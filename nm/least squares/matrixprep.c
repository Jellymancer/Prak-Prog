#include "matrixprep.h"
#include<math.h>


double funs(int i, double x){
   		switch(i){
   			case 0: return log(x); break;
   			case 1: return 1;   break;
   			case 2: return x;     break;
   			default: return NAN;
   			}
		}



void matrix_generation(gsl_matrix* M,gsl_vector* b,double* x,double* y,double* dy){

int n=M->size1;
int m=M->size2;
for(int i=0;i<n;i++){
gsl_vector_set(b,i,*(y+i) / *(dy+i));
for(int j=0;j<m;j++){
	gsl_matrix_set(M,i,j,funs(j,*(x+i)) / *(dy+j));
}}

}


void printm(gsl_matrix* T,FILE* f){
for(int c=0;c< T->size1;c++){
	for(int r=0;r< T->size2;r++){
		fprintf(f,"%7.3f  ",gsl_matrix_get(T,c,r));
}
fprintf(f,"\n");
}
fprintf(f,"\n");
}

void printv(gsl_vector *v,FILE* f){
	for(int r=0;r< v->size;r++){
		fprintf(f,"%7.3f",gsl_vector_get(v,r));
		fprintf(f,"\n");}
}
