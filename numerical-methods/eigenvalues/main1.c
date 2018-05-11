#include "jac.h"

int main(int argc, char** argv){
int n, sweeps;
n=atoi(argv[1]);

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

sweeps = jac(A,e,V);
printf("\n\n\n\nAssigment A1:\n");
printf("Jacobi eigenvalue algorythm implementation\n");
printf("\na Generating a random %dx%d symmetric matrix A\n",n,n); printm(B,stdout);
printf("\nRunning the jacobi diagonalization algorythm gives the following eigenvalues: \n");
printf( "eig=:\n");
printv(e,stdout);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,B,V,0,A);
gsl_blas_dgemm(CblasTrans  ,CblasNoTrans,1,V,A,0,B);
printf( "V^T*A*V should be diagonal with above eigenvalues present in the diagonal:\n V^T*A*V =\n");
printm(B,stdout);
printf("The diagonalization algorythm is run for random matrices of increasing size n. The time taken for the algorythm to run is");
printf(" measured and can be seen in figure 1. It should scale as n³ which is checked by fitting a f(n)=a*n³ to the points.\n");
}
