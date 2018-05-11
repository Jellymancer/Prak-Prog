#include "jac.h"

int main(){
int nmax=14, sweeps;
FILE* swp=fopen("swp.txt","w");
for(int n=2;n<nmax;n++){
gsl_matrix* A=gsl_matrix_calloc(n,n);
gsl_matrix* B=gsl_matrix_calloc(n,n);
gsl_matrix* Acpy=gsl_matrix_calloc(n,n);
srand(21); // seed random generator. Different values produce different matrices.



for(int i=0;i<n;i++){
for(int j=i;j<n;j++){
double x=RND;
gsl_matrix_set(A,i,j,x);
gsl_matrix_set(A,j,i,x);
}}
gsl_matrix_memcpy(B,A);
gsl_matrix_memcpy(Acpy,A);
gsl_matrix *V = gsl_matrix_alloc(n,n);
gsl_vector *e = gsl_vector_alloc(n);

if(n==3){

gsl_matrix_memcpy(A,Acpy);
sweeps = jac(A,e,V);
int row=1; //the desired row to reduce.
int ascending = 1; // specify if the algorythm should go from lowest to highest eigenvalue (1) or vice versa (0)
//Run the algoruthm for thespecified row.
gsl_matrix_memcpy(A,Acpy);
double* result_onevalue = jac_onevalue(A,e,V,row,ascending);

printf("\nAssigment B1:\n");
printf("Jacobi eigenvalue algorythm implementation lowest eigenvalue only.\n");
printf("\nMatrix A from assignment A is used here, which means that the eigenvalues found here should correspond to the ones foudn in the previous assignment.\n");
printf("\nRunning the jacobi diagonalization algorythm for the first row only: \n");
printf( "The lowest eigenvalue is: l=%7.4f:\n",*(result_onevalue));
printf( "\n Assignment B2: The algorythm can be set up to give the highest eigenvalue first, then the second highest etc. by ");
printf( "changing the sign of the rotation angle.\n");
ascending = 0;
gsl_matrix_memcpy(A,Acpy);

double* result_onevalue2 = jac_onevalue(A,e,V,row,ascending);
printf("The highest eigenvalue is: l=%7.4f\n",*(result_onevalue2));
printf("\nAssignment B3: The number of sweeps required for the diagonalization of the entire matrix by the cyclic method is %d compared to %1.0f",sweeps,*(result_onevalue+1));
printf(" sweeps required to find the lowest eigenvalue only\n");
printf("\nAssignment B4: The number of sweeps required to diagonilize increasingly larger matrices using the two methods is shown in figure 2\n");
gsl_matrix_memcpy(A,Acpy);
}

gsl_matrix_memcpy(A,Acpy);
int sweeps_cyc = jac(A,e,V);
gsl_matrix_memcpy(A,Acpy);
double sweeps_eigbyeig=jac_eigvlbyeigvl(A,e,V);
fprintf(swp,"%d %d %g\n",n,sweeps_cyc,sweeps_eigbyeig);
}
fclose(swp);
}


