#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<stdio.h>
#include "print.h"
#include<math.h>
#include "func.h"
// functions and Jacobians are found in functions.c

int main() {
	int ncalls=0;
	int n=2;  //size
	double eps = 1e-7; //tolerance
	gsl_vector* x=gsl_vector_alloc(n);
	gsl_matrix* J=gsl_matrix_alloc(n,n);

	gsl_vector_set(x,0,2);
	gsl_vector_set(x,1,2);
	gsl_vector* fx=gsl_vector_alloc(n);
	printf("Finding the extremum of the Rosenbrock's function:\n");
	printf("initial guess v=[x,y]:\n"); printv(x,stdout);
	ncalls=newton(fRos,Jacobi_fRos,x,1e-3,eps,1);
	printf("The minimum is found at [x,y]:\n"); printv(x,stdout);
	printf("The number of function calls is %i\n",ncalls);
	fRos(x,fx);
	printf("Inserting the found values of x and y into Rosenbrocks' function yields (Should be zero):\n"); printv(fx,stdout);
	gsl_vector_set(x,0,2);
	gsl_vector_set(x,1,2);
	ncalls = newton(fRos,Jacobi_fRos,x,1e-3,eps,0);
	printf("\nUsing the newton algorythm with a maunually supplied analytical jacobian gives a minimum at [x,y]:\n"); printv(x,stdout);
	printf("The number of function calls is %i\n",ncalls);
	fRos(x,fx);
	printf("Inserting the found values of x and y into Rosenbrocks' function yields (Should be zero):\n"); printv(fx,stdout);

//System of equations...
	printf("\n\n");
	ncalls=0;
	gsl_vector_set(x,0,5);
	gsl_vector_set(x,1,20);
	printf("Root finding for the system of equations given in assignment A:\n");
	printf("initial guess v=[x,y]:\n"); printv(x,stdout);
	fSys(x,fx);
	ncalls=newton(fSys,Jacobi_fSys,x,1e-5,eps,1);
	printf("Solution v=[x,y]:\n"); printv(x,stdout);
	printf("The number of function calls is %i\n",ncalls);
	fSys(x,fx);
	printf("Putting the found values of x and y into the system of equations gives (Should be zero):\n"); printv(fx,stdout);

	gsl_vector_set(x,0,5);
	gsl_vector_set(x,1,20);
	ncalls=newton(fSys,Jacobi_fSys,x,1e-5,eps,0);
	printf("\nUsing the newton algorythm using a maunually supplied analytical jacobian gives a minimum at [x,y]:\n");
	printv(x,stdout);
	printf("The number of function calls is %i\n",ncalls);
	fSys(x,fx);
	printf("Putting the found values of x and y into the system of equations gives (Should be zero):\n"); printv(fx,stdout);


//Himmelblau
	printf("\n\n");
	ncalls=0;
	gsl_vector_set(x,0,-32);
	gsl_vector_set(x,1,-52);
	printf("Finding the extremum of Himmelblaus' function:\n");
	printf("initial guess v=[x,y]:\n"); printv(x,stdout);
	ncalls=newton(fHim,Jacobi_fHim,x,1e-6,eps,1);
	printf("The minimum is found at [x,y]:\n"); printv(x,stdout);
	printf("The number of function calls is %i\n",ncalls);
	fHim(x,fx);
	printf("Inserting the found values of x and y into Himmelblaus' function yields (Should be zero):\n"); printv(fx,stdout);

	gsl_vector_set(x,0,-32);
	gsl_vector_set(x,1,-52);
	ncalls=newton(fHim,Jacobi_fHim,x,1e-6,eps,0);
	printf("\nUsing the newton algorythm using a maunually supplied analytical jacobian gives a minimum at [x,y]:\n");printv(x,stdout);
	printf("The number of function calls is %i\n",ncalls);
	fHim(x,fx);
	printf("Inserting the found values of x and y into Himmelblaus' function yields (Should be zero):\n"); printv(fx,stdout);

}
