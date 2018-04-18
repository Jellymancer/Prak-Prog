#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include<math.h>


int main()
{

/* Here I create the matrix and vector from the assignment */
float array[9];
*array=6.13;
*(array+1)=-2.90;
*(array+2)=5.86;
*(array+3)=8.08;
*(array+4)=-6.31;
*(array+5)=-3.89;
*(array+6)=-4.36;
*(array+7)=1.00;
*(array+8)=0.19;



int n=3;
gsl_matrix* A=gsl_matrix_calloc(n,n);

for(int i=0;i<n;i++)
for(int j=0;j<n;j++){
if(i==0){
	gsl_matrix_set(A,i,j,*(array+j));}
else if(i==1){
	gsl_matrix_set(A,i,j,*(array+j+3));}
else{
gsl_matrix_set(A,i,j,*(array+j+6));}}

printf("original matrix M is [M(row,column)]\n");
for(int i=0;i<n;i++)
for(int j=0;j<n;j++)
	printf("m(%d,%d) = %g\n",i,j,gsl_matrix_get(A,i,j));


float q[3];
*q=6.23;
*(q+1)=5.37;
*(q+2)=2.29;

gsl_vector* b=gsl_vector_calloc(n);

for(int i=0;i<n;i++)
	gsl_vector_set(b,i,*(q+i));




/* Now I solve the system of equations*/


gsl_vector* x=gsl_vector_calloc(n);
gsl_linalg_HH_solve(A,b,x);

gsl_matrix* B=gsl_matrix_calloc(n,n);

for(int i=0;i<n;i++)
for(int j=0;j<n;j++){
if(i==0){
	gsl_matrix_set(B,i,j,*(array+j));}
else if(i==1){
	gsl_matrix_set(B,i,j,*(array+j+3));}
else{
gsl_matrix_set(B,i,j,*(array+j+6));}}

gsl_vector* y=gsl_vector_calloc(n);
gsl_blas_dgemv(CblasNoTrans,1,B,x,0,y);

printf("\n \n The solution vector x is found to be [x(index)]:\n");
/* and print */
for( int i=0;i<n;i++)
	printf("x(%d)=%g\n",i,gsl_vector_get(x,i));

printf("\n \n Using x to check if the solution works yields [b(index)]:\n");
for (int i=0;i<n;i++)
	printf("b(%d)= %g\n",i,gsl_vector_get(y,i));



gsl_matrix_free(A);
gsl_matrix_free(B);
return 0;
}
