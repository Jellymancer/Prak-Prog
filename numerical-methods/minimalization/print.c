#include "func.h"

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
