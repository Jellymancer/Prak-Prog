#include "func.h"
void qrgsdecomp(gsl_matrix *E,gsl_matrix *W){
int s=E->size2; //Number of columns in matrix E
for(int i=0;i<s;i++){
        gsl_vector_view col=gsl_matrix_column(E,i);
        double Rii = gsl_blas_dnrm2(&col.vector);
        gsl_matrix_set(W,i,i,Rii);
        gsl_vector_scale(&col.vector,1/Rii);

        for(int j=i+1;j<s;j++){
                gsl_vector_view col2= gsl_matrix_column(E,j);
                double Rij = 0;
                gsl_blas_ddot(&col.vector,&col2.vector,&Rij);
                gsl_blas_daxpy(-Rij,&col.vector,&col2.vector);
                gsl_matrix_set(W,i,j,Rij);
                }}}

void qrgssolve(gsl_matrix* Q,  gsl_matrix* R,gsl_vector* b,gsl_vector* c){
gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,c);
for(int i=c->size-1; i>=0; i--){
        double s=gsl_vector_get(c,i);
        for(int k=i+1;k< c->size; k++)
                s-=gsl_matrix_get(R,i,k)*gsl_vector_get(c,k);
                gsl_vector_set(c,i,s/gsl_matrix_get(R,i,i));}}

int newton(void f(gsl_vector* x,gsl_vector* fx), void Jacobi(gsl_vector* p, gsl_matrix* J),gsl_vector* x, double dx, double eps, int numericJac){
	int n=x->size;
	gsl_matrix* J = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* z  = gsl_vector_alloc(n);
	gsl_vector* fz = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);
	int ncalls = 0;
	do{ ncalls ++; f(x,fx);
		if(numericJac == 1){
		for (int j=0;j<n;j++){
			gsl_vector_set(x,j,gsl_vector_get(x,j)+dx); // advance x by dx
			f(x,df); // df=f(x+dx)
			gsl_vector_sub(df,fx); /* df=f(x+dx)-f(x) */
			for(int i=0;i<n;i++) gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/dx); //for loop for finding numerical jacobian
			gsl_vector_set(x,j,gsl_vector_get(x,j)-dx);
			}}
		else{Jacobi(x,J);} // if analytic Jacobian is supplied, use it instead.

		qrgsdecomp(J,R);
		qrgssolve(J,R,fx,Dx);
		gsl_vector_scale(Dx,-1); // calculate Dx newton step
		double lam=1;
		do{
			gsl_vector_memcpy(z,x);
			gsl_blas_daxpy(lam,Dx,z); // take dx step
			f(z,fz);
			lam/=2.0;
			}while(gsl_blas_dnrm2(fz)>(1-lam/2)*gsl_blas_dnrm2(fx) && lam>0.02); // continue until lambda meets condition
			// || f(x + lam ∆x)|| < (1−lam/2)*||f(x)||
		gsl_vector_memcpy(x,z); // replace previous values with vales after step
		gsl_vector_memcpy(fx,fz);
//		for(int i=0;i<n;i++){assert(isnan(gsl_vector_get(x,i)) != 0 && isinf(gsl_vector_get(x,i)) !=0);}
		}while(gsl_blas_dnrm2(Dx)>dx && gsl_blas_dnrm2(fx)>eps);//Do as long as the x and y steps are under the given tolerances.

	return ncalls;
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(fx);
	gsl_vector_free(z);
	gsl_vector_free(fz);
	gsl_vector_free(df);
	gsl_vector_free(Dx);
}
