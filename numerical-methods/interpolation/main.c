#include "interpol.h"
#include<math.h>

int main(){
	int size = 15;
	double x[size],y[size]; //allocate x and y points and make a sine curve to be interpolated.
	double dydx;
	FILE* sine;
	sine = fopen("sine.dat","w");
	for(int i=0;i<size;i++){
		x[i]=i;
		y[i]=sin(0.5*i);
		dydx=0.5*cos(0.5*i);
		fprintf(sine,"%g %g %g\n",x[i],y[i],dydx);
	}
	fclose(sine);

	FILE* splines;
	splines = fopen("splines.dat","w");
	double dz=0.1; //stepsize
	for(double z=x[0];z<=x[size-1];z+=dz){ // linear spline is calculated.
		double l=linterp(size,x,y,z);
		fprintf(splines,"%g %g\n",z,l);
	}
	fclose(splines);

	double integral;
	double zmax=5; // max value for integration.
	integral = linterpint(size,x,y,zmax);
	printf("\n \n Assignment 1 and 2: Integrations\nFunction f(x)=sin(0.5*x) and the linear and quadratic spline are shown in figure 1\n"); 
	printf("Calculating the and spline integral of the cosine function");
	printf(" from %g to x=%g yields a value of I_lin=%g.\nManual integration of the function gives I_man=3.60229\n",x[0],zmax,integral);


	// Q-spline
	qspline* Q=qspline_alloc(size,x,y);

	FILE* qsplines;
	qsplines=fopen("qspline.dat","w");

	for(double z=x[0];z<=x[size-1];z+=dz){
		double l=qspline_eval(Q,z);
		double dl=qspline_deriv(Q,z);
		fprintf(qsplines,"%g %g %g\n",z,l,dl);
	}
	fclose(qsplines);
	double qintegral = qspline_integral(Q,zmax);
	printf("The same integral performed using the quadratic spline yields I_quad = %g\n",qintegral);
	printf("The integral of the quadratic spline is shown together with the analytic integral in figure 2\n");

	return 0;
}

