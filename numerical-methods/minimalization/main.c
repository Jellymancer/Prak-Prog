#include "func.h"

int main() {
	int ncalls; //no of calls
	int n=2; //size
	double fval;
	gsl_vector* v=gsl_vector_alloc(n);
        gsl_vector* fx=gsl_vector_alloc(n);

	gsl_vector_set(v,0,1.5); //Startværdier
        gsl_vector_set(v,1,3);
	printf("\n\nAssignment A: finding minima of Rosenbrockss' and Himmelblaus' functions using newtons method.\n");
        printf("Finding the minimum of Rosenbrocks' function:\n");
        printf("initial guess v=[x,y]:\n"); printv(v,stdout);
        fRos(v,fx);
        ncalls = newton(fRos,Jacobi_fRos,v,1e-5,1e-7);

        printf("\nThe minimum is found at [x_min,y_min]:\n"); printv(v,stdout);
	fval=fRos(v,fx);
	printf("The function value at this point is f(x_min,y_min)=%g\n",fval);
	printf("\nIt took %d steps to find the minimum\n",ncalls);

	gsl_vector_set(v,0,4); //Startværdier
        gsl_vector_set(v,1,3);
        printf("\n\nFinding one of the minima of Himmelblaus' function:\n");
        printf("initial guess v=[x,y]:\n"); printv(v,stdout);
        fHim(v,fx);
        ncalls = newton(fHim,Jacobi_fHim,v,1e-5,1e-7);
        printf("\nThe minimum is found at [x_min,y_min]:\n"); printv(v,stdout);
	fval=fHim(v,fx);
	printf("The function value at this point is f(x_min,y_min)=%g\n",fval);
	printf("\nIt took %d steps to find the minimum\n",ncalls);

// Broyden update
	printf("\n\nAssignment B: Broyden update.\nTrying the same two function using quasi-newtion with Broyden update.\n");
	gsl_vector_set(v,0,1.5);
        gsl_vector_set(v,1,3);
        printf("Finding the minimum of Rosenbrocks' function using broyden:\n");
        printf("initial guess v=[x,y]:\n"); printv(v,stdout);
        fRos(v,fx);
        ncalls = newton_broy(fRos,v,1e-5);
        printf("\nThe minimum is found at [x,y]:\n"); printv(v,stdout);
	fval=fRos(v,fx);
	printf("The function value at this point is f(x_min,y_min)=%g\n",fval);
	printf("\nIt took %d steps to find the minimum\n",ncalls);

	gsl_vector_set(v,0,3.5); //Startværdier
        gsl_vector_set(v,1,1);
        printf("\n\nFinding one of the minima of Himmelblaus' function using broyden update:\n");
        printf("initial guess v=[x,y]:\n"); printv(v,stdout);
        fHim(v,fx);
        ncalls = newton_broy(fHim,v,1e-5);
        printf("\nThe minimum is found at [x_min,y_min]:\n"); printv(v,stdout);
	fval=fHim(v,fx);
	printf("The function value at this point is f(x_min,y_min)=%g\n",fval);
	printf("\nIt took %d steps to find the minimum\n",ncalls);


	gsl_vector_free(v); gsl_vector_free(fx);
// Exponential decay. Data listed in exp.dat.


/*	double t[]= {0.47,1.41,2.36,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
	double y[]= {5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
	double e[]= {0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};
	int m = sizeof(t)/sizeof(t[0]);
*/ // For some reason this dosent converge at all. using a modified data set ...

double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int m = sizeof(t)/sizeof(t[0]);

	printf("\n\nRadioactive decay.\nNow I make a fit to the radioactive data (listed in exp.dat) using least squares.\n");
	n=3;
	gsl_vector* v2=gsl_vector_alloc(n);
	gsl_vector_set(v2,0,5.13); //Startværdier
        gsl_vector_set(v2,1,2.7);
        gsl_vector_set(v2,2,1.08);

        ncalls = newton_broy(fRad,v2,1e-5);
        printf("\nThe parameters of the best fit are found to be: [A,T,B]:\n"); printv(v2,stdout);
	printf("The data points together with the best fit can be seen in figure 1.\n");


double Amin=gsl_vector_get(v2,0), Tmin=gsl_vector_get(v2,1),Bmin=gsl_vector_get(v2,2);
FILE *fp;
fp = fopen("exp.dat","w");
for(int i=0;i<10;i++) fprintf(fp,"%g\t%g\t%g\n", t[i],y[i],e[i]);
fprintf(fp,"\n\n");
for(double j=0.0; j<=9;j+=0.1) fprintf(fp, "%g %g\n", j, Amin*exp(-(j)/Tmin)+Bmin);

fclose(fp);

	gsl_vector_free(v2);
return 0;
}
