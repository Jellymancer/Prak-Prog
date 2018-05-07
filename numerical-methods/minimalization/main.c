#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<stdio.h>
#include "print.h"
#include<math.h>
#include "func.h"

int main() {
	int ncalls; //no of calls
	int n=2; //size

	gsl_vector* v=gsl_vector_alloc(n);
	gsl_vector* df=gsl_vector_alloc(n);
	gsl_matrix* H=gsl_matrix_alloc(n,n);

	gsl_vector_set(v,0,-2.4); //Startværdier
        gsl_vector_set(v,1,0.4);

        gsl_vector* fx=gsl_vector_alloc(n);
        printf("Finding the extremum of Rosenbrocks' function:\n");
        printf("initial guess v=[x,y]:\n"); printv(v,stdout);
        fRos(v,fx);

        ncalls = newton(fRos,Jacobi_fRos,v,1e-5,1e-7,0);
        printf("\nThe minimum is found at [x,y]:\n"); printv(v,stdout);
	printf("\nIt took %d steps to find the minimum\n",ncalls);

	gsl_vector_set(v,0,-2); //Startværdier
        gsl_vector_set(v,1,3);
        printf("\n\nFinding the extremum of Himmelblaus' function:\n");
        printf("initial guess v=[x,y]:\n"); printv(v,stdout);
        fHim(v,fx);
        ncalls = newton(fHim,Jacobi_fHim,v,1e-5,1e-7,0);
        printf("\nThe minimum is found at [x,y]:\n"); printv(v,stdout);
	printf("\nIt took %d steps to find the minimum\n",ncalls);


return 0;
}
