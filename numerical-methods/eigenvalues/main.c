#include<stdio.h>
#include<stdlib.h>
#define RND (double)rand()/RAND_MAX
#include<time.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


void printm(gsl_matrix*,FILE*);
void printv(gsl_vector*,FILE*);
int jac(gsl_matrix* , gsl_vector*,gsl_matrix*);
double* jac_onevalue(gsl_matrix* , gsl_vector*,gsl_matrix*,int,int);
double jac_eigvlbyeigvl(gsl_matrix* , gsl_vector*,gsl_matrix*);
double jac_classic(gsl_matrix*, gsl_vector*, gsl_matrix*);

int main(){
int nmax=10;
for(int n=1;n<nmax;n++){

gsl_matrix* A=gsl_matrix_calloc(n,n);
gsl_matrix* B=gsl_matrix_calloc(n,n);
gsl_matrix* Acpy=gsl_matrix_calloc(n,n);

srand(21);



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

int sweeps = jac(A,e,V);
if(n==3){
FILE* f;
f=fopen("opgA.txt","w");
fprintf(f,"Assigment A1:\n");
fprintf(f,"Jacobi eigenvalue algorythm implementation\n");
fprintf(f,"\na Generating a random %dx%d symmetric matrix A\n",n,n); printm(B,f);
fprintf(f,"\nRunning the jacobi diagonalization algorythm gives the following eigenvalues: \n");
fprintf(f, "eigenvalues:\n");
printv(e,f);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,B,V,0,A);
gsl_blas_dgemm(CblasTrans  ,CblasNoTrans,1,V,A,0,B);
fprintf(f, "check: V^T*A*V should be diagonal with above eigenvalues:\n");
printm(B,f);
fclose(f);

gsl_matrix_memcpy(A,Acpy);
FILE* f2;
f2=fopen("opgB.txt","w");
int row=1; //the desired row to reduce.
int ascending = 1; // specify if the algorythm should go from lowest to highest eigenvalue (1) or vice versa (0)
//Run the algoruthm for thespecified row.
double* result_onevalue = jac_onevalue(A,e,V,row,ascending);

fprintf(f2,"Assigment B1:\n");
fprintf(f2,"Jacobi eigenvalue algorythm implementation lowest eigenvalue only.\n");
fprintf(f2,"\nMatrix A from assignment A is used here, which means that the eigenvalues found here should correspond to the ones foudn in the previous assignment.\n");
fprintf(f2,"\nRunning the jacobi diagonalization algorythm for the first row only: \n");
fprintf(f2, "The lowest eigenvalue is: l=%7.4f:\n",*(result_onevalue));
fprintf(f2, "\n Assignment B2: The algorythm can be set up to give the highest eigenvalue first, then the second highest etc. by ");
fprintf(f2, "changing the sign of the rotation angle.\n");
ascending = 0;
gsl_matrix_memcpy(A,Acpy);

double* result_onevalue2 = jac_onevalue(A,e,V,row,ascending);
fprintf(f2,"The highest eigenvalue is: l=%7.4f\n",*(result_onevalue2));
fprintf(f2,"\nAssignment B3: The number of sweeps required for the diagonalization of the entire matrix by the cyclic method is %d compared to %1.0f",sweeps,*(result_onevalue+1));
fprintf(f2,"sweeps required to find the lowest eigenvalue only\n");
fprintf(f2,"\nAssignment B4: The number of sweeps required to diagonilize the matix eigenvalue using the two methods is shown in the attached figure\n");
gsl_matrix_memcpy(A,Acpy);
fclose(f2);
}

gsl_matrix_memcpy(A,Acpy);
double totalsweeps=jac_eigvlbyeigvl(A,e,V);
gsl_matrix_memcpy(A,Acpy);
double classicsweeps = jac_classic(A,e,V);
printf("%d %d %1.0f %1.0f\n",n,sweeps,totalsweeps,classicsweeps);
}
}


