#include <stdlib.h>
#include <assert.h>
#include "interpol.h"


cspline* cspline_alloc(int n,double* x,double* y){ //builds cspline
	cspline* s = malloc(sizeof(cspline));//spline is generated. 
	s->x = (double*) malloc(n*sizeof(double));
	s->y = (double*) malloc(n*sizeof(double));
	s->b = (double*) malloc(n*sizeof(double));
	s->c = (double*) malloc((n-1)*sizeof(double));
	s->d = (double*) malloc((n-1)*sizeof(double));
	s->n = n;

	for(int i=0;i<n;i++){//x and y lists are defined in the c-spline
		s->x[i]=x[i];
		s->y[i]=y[i];}

	double h[s->n]; double dy[s->n]; double p[s->n];
	for (int i=0;i<n-1;i++){
		h[i]=x[i+1]-x[i];
		dy[i]=y[i+1]-y[i];
		p[i]=dy[i]/h[i];}

	//Tridiagonal system is set up.
	double D[s->n], Q[s->n-1], B[s->n];
	D[0]=2, Q[0]=1;
	for (int i=0;i<n-2;i++){
		D[i+1]=2*h[i]/h[i+1]+2;
		D[n-1]=2;
		Q[i+1]=h[i]/h[i+1];
		B[i+1]=3*(p[i]+p[i+1]*h[i]/h[i+1]);}

	//Now we use gauss elimination of the system
	B[0]=3*p[0];
	B[n-1]= 3*p[n-2];
	for (int i=1;i<n;i++){
		D[i]-=Q[i-1]/D[i-1];
		B[i]-=B[i-1]/D[i-1];}

	//Back substitution
	s->b[n-1]=B[n-1]/D[n-1];
	for (int i=n-2;i>=0;i--){
		s->b[i]=(B[i]-Q[i]*s->b[i+1])/D[i];}
	for (int i=0;i<n-1;i++){
		s->c[i]=(-2*s->b[i]-s->b[i+1]+3*p[i])/h[i];
		s->d[i]=(s->b[i]+s->b[i+1]-2*p[i])/(h[i]*h[i]);}
	return s;

}
double cspline_eval(cspline * s, double z){
	int i = binsearchcspl(s->n,z,s);
	double h=z-s->x[i];
	return s->y[i]+h*(s->b[i]+h*(s->c[i]+h*s->d[i]));
}

double cspline_deriv(cspline* s, double z){//calculates the derivative of the cspline.
	int i = binsearchcspl(s->n,z,s);
	double h=z-s->x[i];
	return s->b[i]+2*s->c[i]*h+3*s->d[i]*h*h;
}

double cspline_deriv2(cspline* s, double z){//calculates the second derivative of the cspline.
	int i = binsearchcspl(s->n,z,s);
	double h=z-s->x[i];
	return 2*s->c[i]+6*s->d[i]*h;
}

double cspline_integ(cspline* s, double z){//calculates the derivative of the cspline.
	int i = binsearchcspl(s->n,z,s);

	double sum=0;
	for(int k=0;k<i;k++){//integral for each interval is found and summed
		double xk=s->x[k], xk1=s->x[k+1], bk=s->b[k], ck=s->c[k];
		double yk=s->y[k], dk=s->d[k], dx=xk1-xk;
		sum=sum+yk*dx+bk*0.5*dx*dx+ ck*dx*dx*dx/3+ dk*dx*dx*dx*dx/4;}
	double dxi=z-s->x[i];
	sum+= s->y[i]*dxi+s->b[i]*0.5*dxi*dxi + s->c[i]*dxi*dxi*dxi/3 + s->d[i]*dxi*dxi*dxi*dxi/4;
	return sum;
}


void cspline_free(cspline * s){
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s->d);
	free(s);}
