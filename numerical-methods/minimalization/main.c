#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<stdio.h>
#include "print.h"
#include<math.h>

int newton(
	void f(gsl_vector* x,gsl_vector* fx,gsl_matrix* J),
	gsl_vector* x, double dx, double eps,gsl_matrix* J);

int main() {
	int ncalls=0;
	int n=2;
	gsl_vector* v=gsl_vector_alloc(n);
	gsl_vector* df=gsl_vector_alloc(n);
	gsl_matrix* H=gsl_matrix_alloc(n,n);

	void fRos(gsl_vector* v,gsl_vector* df, gsl_matrix* H){
		ncalls++;
		double x=gsl_vector_get(v,0), y=gsl_vector_get(v,1);
		gsl_vector_set(df,0, 2*(1-x)*(-1)+100*2*(y-x*x)*(-1)*2*x);
		gsl_vector_set(df,1, 100*2*(y-x*x));
		gsl_matrix_set(H,0,0,2+100*2*(-2*x)*(-2*x)+100*2*(y-x*x)*(-1)*2); gsl_matrix_set(H,0,1,100*2*(-2*x));
		gsl_matrix_set(H,1,0,800*x*y); gsl_matrix_set(H,1,1,100*2);
		}

	void fHim(gsl_vector* v, gsl_vector* df, gsl_matrix* H){
		ncalls++;
		double x=gsl_vector_get(v,0), y=gsl_vector_get(v,1);
                gsl_vector_set(df,0,4*x*x*x+2*x*(2*y-21)+2*(y*y-7));
		gsl_vector_set(df,1,4*y*y*y+2*y*(2*x-13)+2*x*x-22);
		gsl_matrix_set(H,0,0,4*3*x*x+2*(2*y-21)); gsl_matrix_set(H,0,1,2*y*2+2*2*x);
		gsl_matrix_set(H,1,0,2*x*2+2*2*y); gsl_matrix_set(H,1,1,4*3*y*y+2*(2*x-13));
		}

	gsl_vector_set(v,0,-2); //Startværdier
        gsl_vector_set(v,1,8);
        gsl_vector* fx=gsl_vector_alloc(n);
        printf("Finding the extremum of the Rosenbrock's function:\n");
        printf("initial guess v=[x,y]:\n"); printv(v,stdout);
        fRos(v,fx,H);
        printf("\nAnd the value of the derivatives at the start Df=[dfdx,dfdy]:\n");
	printv(fx,stdout);
        newton(fRos,v,1e-6,1e-3,H);
        printf("\nThe minimum is found at [x,y]:\n"); printv(v,stdout);


return 0;
}
