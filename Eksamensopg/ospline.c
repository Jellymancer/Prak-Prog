#include <stdlib.h>
#include <assert.h>
#include "interpol.h"


ospline* ospline_alloc(int n,double* x,double* y, double* df){ //builds qspline
	ospline* s = malloc(sizeof(ospline));//spline is generated. 
	s->x = (double*) malloc(n*sizeof(double));
	s->y = (double*) malloc(n*sizeof(double));
	s->df = (double*) malloc(n*sizeof(double));
	s->c = (double*) malloc((n-1)*sizeof(double));
	s->d = (double*) malloc((n-1)*sizeof(double));
	s->n = n;

	for(int i=0;i<n;i++){//x,y and df lists are defined in the q-spline (df being the derivative)
		s->x[i]=x[i];
		s->y[i]=y[i];
		s->df[i]=df[i];}

	double h[s->n]; double dy[s->n]; double p[s->n];
	for (int i=0;i<n-1;i++){
		h[i]=x[i+1]-x[i];
		dy[i]=y[i+1]-y[i];}


	for(int i=0;i<n-1;i++)
	s->d[i]=pow(1+1/3,-1)*( (df[i+1]-df[i])*pow(3*pow(h[i],2),-1) - pow(3*pow(h[i],3),-1)*(dy[i]-df[i]*h[i]));

        for(int i=0;i<n-1;i++)  s->c[i]=(dy[i]-df[i]*h[i] - s->d[i]*pow(h[i],3))*pow(h[i],-2);

	return s;

}
double ospline_eval(ospline * s, double z){
	int i = binsearchospl(s->n,z,s);
	double h=z-s->x[i];
	return s->y[i]+h*s->df[i]+h*h*s->c[i]+h*h*h*s->d[i];
}

double ospline_deriv(ospline* s, double z){//calculates the derivative of the ospline.
	int i = binsearchospl(s->n,z,s);
	double h=z-s->x[i];
	return s->df[i]+2*s->c[i]*h+3*s->d[i]*h*h;
}

double ospline_deriv2(ospline* s, double z){//calculates the second derivative of the ospline.
	int i = binsearchospl(s->n,z,s);
	double h=z-s->x[i];
	return 2*s->c[i]+6*s->d[i]*h;
}

double ospline_integ(ospline* s, double z){//calculates the derivative of the cspline.
	int i = binsearchospl(s->n,z,s);

	double sum=0;
	for(int k=0;k<i;k++){//integral for each interval is found and summed
		double xk=s->x[k], xk1=s->x[k+1], bk=s->df[k], ck=s->c[k];
		double yk=s->y[k], dk=s->d[k], dx=xk1-xk;
		sum=sum+yk*dx+bk*0.5*dx*dx+ ck*dx*dx*dx/3+ dk*dx*dx*dx*dx/4;}
	double dxi=z-s->x[i];
	sum+= s->y[i]*dxi+s->df[i]*0.5*dxi*dxi + s->c[i]*dxi*dxi*dxi/3 + s->d[i]*dxi*dxi*dxi*dxi/4;
	return sum;
}


void ospline_free(ospline * s){
	free(s->x);
	free(s->y);
	free(s->df);
	free(s->c);
	free(s->d);
	free(s);}
