#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<stdio.h>
#include "print.h"
#include<math.h>

int newton (
	void f(gsl_vector* x,gsl_vector* fx,gsl_matrix* J),
	gsl_vector* x, double dx, double eps);
int newtonana (
	void f(gsl_vector* x,gsl_vector* fx,gsl_matrix* J),
	gsl_vector* x, double dx, double eps,gsl_matrix* J);

int main() {
	int ncalls=0;
	int n=2;
	gsl_vector* x=gsl_vector_alloc(n);
	gsl_matrix* J=gsl_matrix_alloc(n,n);

	void fRos(gsl_vector* p,gsl_vector* fx, gsl_matrix* J){
		ncalls++;
		double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1);
		gsl_vector_set(fx,0, 2*(1-x)*(-1)+100*2*(y-x*x)*(-1)*2*x);
		gsl_vector_set(fx,1, 100*2*(y-x*x));
		gsl_matrix_set(J,0,0,2+100*2*(-2*x)*(-2*x)+100*2*(y-x*x)*(-1)*2); gsl_matrix_set(J,0,1,100*2*(-2*x));
		gsl_matrix_set(J,1,0,800*x*y); gsl_matrix_set(J,1,1,100*2);
		}
	void fSys(gsl_vector* p, gsl_vector* fx, gsl_matrix* J){
		ncalls++;
		double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1); double A=10000;
                gsl_vector_set(fx,0,1+1/A-exp(-x)-exp(-y));
		gsl_vector_set(fx,1,A*x*y-1);
		gsl_matrix_set(J,0,0,exp(-x)); gsl_matrix_set(J,0,1,A*y); gsl_matrix_set(J,1,0,exp(-y)); gsl_matrix_set(J,1,1,A*x);
		}

	void fHim(gsl_vector* p, gsl_vector* fx, gsl_matrix* J){
		ncalls++;
		double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1);
                gsl_vector_set(fx,0,4*x*x*x+2*x*(2*y-21)+2*(y*y-7));
		gsl_vector_set(fx,1,4*y*y*y+2*y*(2*x-13)+2*x*x-22);
		gsl_matrix_set(J,0,0,4*3*x*x+2*(2*y-21)); gsl_matrix_set(J,0,1,2*y*2+2*2*x);
		gsl_matrix_set(J,1,0,2*x*2+2*2*y); gsl_matrix_set(J,1,1,4*3*y*y+2*(2*x-13));
		}
	gsl_vector_set(x,0,-2);
	gsl_vector_set(x,1,8);
	gsl_vector* fx=gsl_vector_alloc(n);
	printf("Finding the extremum of the Rosenbrock's function:\n");
	printf("initial guess v=[x,y]:\n"); printv(x,stdout);
	fRos(x,fx,J);
	printf("\nAnd the value of the derivatives at the start Df=[dfdx,dfdy]:\n"); printv(fx,stdout);
	newton(fRos,x,1e-6,1e-3);
	printf("\nThe minimum is found at [x,y]:\n"); printv(x,stdout);
	fRos(x,fx,J);
	newtonana(fRos,x,1e-8,1e-5,J);
	printf("\nUsing the newton algorythm using a maunually supplied analytical jacobian gives a minimum at [x,y]:\n"); printv(x,stdout);

//System of equations...
	printf("\n\n");
	ncalls=0;
	gsl_vector_set(x,0,0.134);
	gsl_vector_set(x,1,10);
	printf("Root finding for the system of equation in opg A:\n");
	printf("initial guess v=[x,y]:\n"); printv(x,stdout);
	fSys(x,fx,J);
	newton(fSys,x,1e-8,1e-5);
	printf("\nsolution v=[x,y]:\n"); printv(x,stdout);
	fSys(x,fx,J);
	printf("\nThe value of the equations using the found x and y. (Should be zero):\n"); printv(fx,stdout);
	fSys(x,fx,J);
	newtonana(fSys,x,1e-8,1e-5,J);
	printf("\nUsing the newton algorythm using a maunually supplied analytical jacobian gives a minimum at [x,y]:\n"); printv(x,stdout);

//Himmelblau
	printf("\n\n");
	gsl_vector_set(x,0,-3);
	gsl_vector_set(x,1,-3);
	printf("Finding the extremum of the Himmelblaus' function:\n");
	printf("initial guess v=[x,y]:\n"); printv(x,stdout);
	fHim(x,fx,J);
	printf("\nAnd the value of the derivatives at the start Df=[dfdx,dfdy]:\n"); printv(fx,stdout);
	newton(fHim,x,1e-6,1e-3);
	printf("\nThe minimum is found at [x,y]:\n"); printv(x,stdout);
	fHim(x,fx,J);
	newtonana(fHim,x,1e-8,1e-5,J);
	printf("\nUsing the newton algorythm using a maunually supplied analytical jacobian gives a minimum at [x,y]:\n"); printv(x,stdout);

}
