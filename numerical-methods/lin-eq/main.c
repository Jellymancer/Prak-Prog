#include<stdio.h>
#include<stdlib.h>
#define RND (double)rand()/RAND_MAX
#include<time.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void qr_gs_decomp(gsl_matrix*, gsl_matrix*);
void qrgssolve(gsl_matrix*,gsl_matrix*,gsl_vector*,gsl_vector*);
void printm(gsl_matrix*);
void printv(gsl_vector*);

int main(){

int n=4;// Definition of nxm matrixes
int m=3;
gsl_matrix* M=gsl_matrix_calloc(n,m);
gsl_matrix* R=gsl_matrix_calloc(m,m);
gsl_matrix* C=gsl_matrix_calloc(n,m);//Matrixes for checking the results.
gsl_matrix* Q=gsl_matrix_calloc(n,m);
gsl_matrix* Qt=gsl_matrix_calloc(m,n);
gsl_matrix* qq=gsl_matrix_calloc(m,m);
srand(time(NULL));

for(int i=0;i<n;i++){
for(int j=0;j<m;j++){
gsl_matrix_set(M,i,j,RND);
}}
printf("Assigment A1:\n");
printf("QR-decomposition of matrix A into matrices Q and R\n");

printf("Random matrix A:\n");
printm(M);
printf("\n\n");


gsl_matrix_memcpy(Q,M);
qr_gs_decomp(Q,R);


printf("QR-decomposition is performed yielding matrix Q: \n");
printm(Q);
printf("\n\n");
printf("and matrix R:\n");
printm(R);

gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Q, R,0.0, C);
gsl_matrix_transpose_memcpy(Qt,Q);
gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Qt, Q,0.0, qq);

printf("\n\nChecking if QR=A.\n QR=\n");
printm(C);

printf("\n\nChecking of QtQ=I. Here Qt is the traspose of Q and I is the identity matrix.\n QtQ=\n");
printm(qq);

// Now to opgA2............................................
gsl_vector* x=gsl_vector_alloc(m);
gsl_vector* bcheck=gsl_vector_alloc(m);
gsl_vector* b=gsl_vector_alloc(m);
gsl_matrix* Mm=gsl_matrix_calloc(m,m);
gsl_matrix* Rr=gsl_matrix_calloc(m,m);
gsl_matrix* Qq=gsl_matrix_calloc(m,m);
gsl_matrix* Qqt=gsl_matrix_calloc(m,m);
gsl_matrix* Mmi=gsl_matrix_calloc(m,m);

for(int i=0;i<m;i++){
for(int j=0;j<m;j++){
gsl_matrix_set(Mm,i,j,RND);
gsl_vector_set(b,i,RND);
}}

printf("\n\n Assignment A2................................................\n");
printf("Randomly generated vector b.\n b=\n");
printv(b);
printf("\n and matrix M.\n M=\n");
printm(Mm);

gsl_matrix_memcpy(Qq,Mm);
qr_gs_decomp(Qq,Rr);
gsl_matrix_transpose_memcpy(Qqt,Qq);
qrgssolve(Qqt,Rr,b,x);

printf("\nNow QR decomposition is performed using the method showcased in A1\n");
printf("Solving R*x=Qt*b using back substitution yields the x-vector.\n x=\n");
printv(x);

printf("\nChecking if Mx equals b.\n Mx=\n");
gsl_blas_dgemv(CblasNoTrans, 1.0, Mm,x, 0.0, bcheck);
printv(bcheck);


//now to Assignment B:
printf("\nAssignment B:.................................\n");
printf("Using the same random matrix M as in A2 I now find its inverse.\nThis is done by using the methods from A1 and A2");
printf(" to solve Mx_i=e_i,\ne_i being the i-th unit vector and x_i being the i-th column of the inverse matrix M⁻¹.\n");

gsl_vector_set_zero(b);
gsl_vector_set_zero(x);
for(int i=0;i< Mm->size2;i++){
	gsl_vector_set(b,i,1);
	qrgssolve(Qqt,Rr,b,x);
	gsl_matrix_set_col(Mmi,i,x);
	gsl_vector_set_zero(b);
	gsl_vector_set_zero(x);
}

printf("The found inverse matrix M⁻¹: \nM⁻¹=\n");
printm(Mmi);
printf("\nChecking if MM⁻¹ = I.\nMM⁻¹=\n");
gsl_matrix_set_zero(Rr);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Mm,Mmi,0.0,Rr);
printm(Rr);

printf("\n\n");





gsl_matrix_free(M);
gsl_matrix_free(R);
gsl_matrix_free(C);
gsl_matrix_free(Qt);
gsl_matrix_free(qq);
gsl_matrix_free(Mm);
gsl_matrix_free(Rr);
gsl_matrix_free(Qqt);
gsl_vector_free(x);
gsl_vector_free(b);
gsl_vector_free(bcheck);
gsl_matrix_free(Mmi);
return 0;
}
