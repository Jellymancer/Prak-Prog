#include <stdlib.h>
#include <assert.h>
#include "interpol.h"


qspline* qspline_alloc(int n,double* x,double* y){ //builds qspline
	qspline* s = malloc(sizeof(qspline));//spline
	s->b = malloc((n-1)*sizeof(double));  // b_i
	s->c = malloc((n-1)*sizeof(double));  // c_i
	s->x = malloc(n*sizeof(double));      // x_i
	s->y = malloc(n*sizeof(double));      // y_i
	s->n = n;
	for(int i=0;i<n;i++){
		s->x[i]=x[i];
		s->y[i]=y[i];
	}
	int i;
	double p[n-1], h[n-1];
	for(i=0;i<n-1;i++){
		h[i]=x[i+1]-x[i];
		p[i]=(y[i+1]-y[i])/h[i];
	}
	s->c[0]=0;

	for(i=0;i<n-2;i++)
		s->c[i+1]=(p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];
	s->c[n-2]/=2;                                 //recursion down:
	for(i=n-3;i>=0;i--)
		s->c[i]=(p[i+1]-p[i]-s->c[i+1]*h[i+1])/h[i];
	for(i=0;i<n-1;i++)
		s->b[i]=p[i]-s->c[i]*h[i];
	return s;
}

double qspline_eval(qspline *s, double z){     //evaluates s(z)
	int n1=s->n, i;
	i=binsearchqspl(n1,z,s);

	double h=z-s->x[i];
	return s->y[i]+h*(s->b[i]+h*s->c[i]);
}

void qspline_free(qspline *s){ //free the allocated memory
	free(s->x); free(s->y); free(s->b); free(s->c); free(s);
}

double qspline_deriv(qspline *s, double z){
	int n1=s->n, i;
	i=binsearchqspl(n1,z,s);

  return s->b[i] +2*s->c[i]*(z-s->x[i]);
}

double qspline_integral(qspline *q, double z){
	int n1=q->n, i;
	i=binsearchqspl(n1,z,q);

  double sum = 0;
  for (int k = 0; k < i; k++) {
	double xk=q->x[k], xk1=q->x[k+1], bk= q->b[k], ck=q->c[k], yk=q->y[k], yk1=q->y[k+1];
	double dx = xk1-xk;
    sum = sum + yk*dx + 0.5*bk*dx*dx + ck*dx*dx*dx/3.0;
  }
sum += q->y[i]*(z-q->x[i]) + 0.5*q->b[i]*(z-q->x[i])*(z-q->x[i]) + q->c[i]*(z-q->x[i])*(z-q->x[i])*(z-q->x[i])/3.0;
  return sum;
}
