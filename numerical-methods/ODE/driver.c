#include "functions.h"
#include<stdlib.h>

int driver(
        double x,                             /* the current value of the variable */
        double x_end,                              /* the end-point of the integration */
        double step,                             /* the current step-size */
        gsl_vector* yx,                             /* the current y(t) */
        double abstol,                            /* absolute accuracy goal */
        double reltol,                            /* relative accuracy goal */
        int max,
	void rkstepX(                          /* the stepper function to be used */
                double x, double h, gsl_vector* yx,
                void f(double x,gsl_vector* yx,gsl_vector* dydx),
                gsl_vector* yxh, gsl_vector* interr),
        void f(double x,gsl_vector *yx,gsl_vector *dydx) /* right-hand-side */
){

int n=yx->size;
gsl_vector* yxh=gsl_vector_alloc(n);
gsl_vector* interr=gsl_vector_alloc(n);


double x_start=x;
int iter=0; // no of iterations.
while(x<x_end){
	if(x+step>x_end) step=x_end-x; // makes the last step fit the end point.

        rkstepX(x,step,yx,f,yxh,interr);

	double err_i=gsl_blas_dnrm2(interr);
	double tau_i=gsl_blas_dnrm2(yxh);

        double tol=(tau_i*reltol+abstol)*sqrt(step/(x_end-x_start)); // local tolerance

        if(err_i<tol){
                iter++;
                if(iter>max-1){printf("ERROR:max iterations reached\n"); break;}
		x=x+step;
	        gsl_vector_memcpy(yx,yxh);}

        if(err_i>0) step*=pow(tol/err_i,0.25)*0.95; else step*=2;

}
gsl_vector_free(yxh);
gsl_vector_free(interr);
return iter+1;
}


int driverwpath(
        double x,                             /* the current value of the variable */
        double x_end,                              /* the end-point of the integration */
        double step,                             /* the current step-size */
        gsl_vector* yx,                             /* the current y(t) */
        double abstol,                            /* absolute accuracy goal */
        double reltol,                            /* relative accuracy goal */
	int max,
        void rkstepX(                          /* the stepper function to be used */
                double x, double h, gsl_vector* yx,
                void f(double x,gsl_vector* yx,gsl_vector* dydx),
                gsl_vector* yxh, gsl_vector* interr),
        void f(double x,gsl_vector *yx,gsl_vector *dydx), /* right-hand-side */
	gsl_matrix* xypath  /* Storage for x and y values */
){

int n=yx->size;
gsl_vector* yxh=gsl_vector_alloc(n);
gsl_vector* interr=gsl_vector_alloc(n);

double x_start=x, err_i, tau_i;
int iter=0; // no of iterations.
while(x < x_end){
	if(x+step>x_end) step=x_end-x; // makes the last step fit the end point.

        rkstepX(x,step,yx,f,yxh,interr);

	err_i=gsl_blas_dnrm2(interr);
	tau_i=gsl_blas_dnrm2(yxh);

        double tol=(tau_i*reltol+abstol)*sqrt(step/(x_end-x_start)); // local tolerance

        if(err_i<tol){
		x=x+step;
        	gsl_vector_memcpy(yx,yxh);
		gsl_matrix_set(xypath,iter,0,x);
		for(int i=1; i<n+1;i++) gsl_matrix_set(xypath,iter,i,gsl_vector_get(yx,i-1));
                iter++;
                if(iter>max-1){printf("ERROR:max iterations reached\n"); break ;}

	}

	if(err_i>0) step*=pow(tol/err_i,0.25)*0.95; else step*=2;

}
gsl_vector_free(yxh);
gsl_vector_free(interr);
return iter;
}

