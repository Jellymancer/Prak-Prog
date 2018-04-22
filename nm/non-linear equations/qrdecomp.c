#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


void qr_gs_decomp(gsl_matrix *E,gsl_matrix *W){

int s=E->size2; //Number of columns in matrix E
for(int i=0;i<s;i++){
	gsl_vector_view col=gsl_matrix_column(E,i);
	double Rii = gsl_blas_dnrm2(&col.vector); // the &-notation is used in order to assure
//that the vector view dosen't get out of scope.

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

