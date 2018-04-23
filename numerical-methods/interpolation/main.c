#include "interpol.h"
#include<math.h>

int main(){
	int size = 15;
	double x[size],y[size];
	FILE* random;
	random = fopen("random.dat","w");
	for(int i=0;i<size;i++){
		x[i]=0.5*i;
		y[i]=sin(0.5*i);
		fprintf(random,"%g %g\n",x[i],y[i]);
	}
	fclose(random);



	FILE* splines;
	splines = fopen("splines.dat","w");

	double dz=0.1;
	for(double z=x[0];z<=x[size-1];z+=dz){
		double l=linterp(size,x,y,z);
		fprintf(splines,"%g %g\n",z,l);
	}
	fclose(splines);

	double integral;
	double zmax=5;
	integral = linterpint(size,x,y,zmax);
	printf("\n \n Assignment 1 and 2: Integrations\n Calculating the integral of the linear spline of the cosine function");
	printf(" from %g to x=%g yields a value of I_lin=-%g.\nManual integration of the function gives I_man=-0.958924\n",x[0],zmax,integral);


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


	return 0;
}

