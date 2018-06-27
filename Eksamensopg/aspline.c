#include <stdlib.h>
#include <assert.h>
#include "interpol.h"

//formula references F(num) are to the interpolation.pdf book from the course.

aspline* aspline_alloc(int n,double* x,double* y){//builds spline
	aspline* s = malloc(sizeof(aspline));//allocating akima spline
	assert(n>2);
	double h[n-1],p[n-1];

	for(int i=0;i<n-1;i++){//set h_i p_i, h_i being x_i+1-x_i.
		h[i]=x[i+1]-x[i];
		assert(h[i]>0);
		p[i]=(y[i+1]-y[i])/h[i];}
//allocating coefficients
	s->Ai = malloc(n*sizeof(double));  // A_i
	s->c = malloc((n-1)*sizeof(double));  // c_i
	s->d = malloc((n-1)*sizeof(double));  // d_i
	s->x = malloc(n*sizeof(double));      // x_i
	s->y = malloc(n*sizeof(double));      // y_i
	s->n = n;

	for(int i=0;i<n;i++){ // set x,y list values in s.
		s->x[i]=x[i];
		s->y[i]=y[i];}
	//two first and last points are determined seperately F(34).
	s->Ai[0]=p[0];  s->Ai[1]=(p[0]+p[1])/2;
	s->Ai[n-1]=p[n-2]; s->Ai[n-2]=(p[n-2]+p[n-3])/2;

	for(int i=2;i<n-2;i++){//finding weigths. F(33)
		double w1=fabs(p[i+1]-p[i]), w2=fabs(p[i-1]-p[i-2]);
		if(w1+w2==0) s->Ai[i]=(p[i-1]+p[i])/2; //determine which formula for A_i' to use F(31/32).
		else s->Ai[i]=(w1*p[i-1]+w2*p[i])/(w1+w2);}

	// Determine c_i and d_i F(29/30).
	for(int i=0;i<n-1;i++){
		s->c[i]=(3*p[i]-2*s->Ai[i]-s->Ai[i+1])/h[i];
		s->d[i]=(s->Ai[i+1]+s->Ai[i]-2*p[i])/h[i]/h[i];
	}
	return s;
}

double aspline_eval(aspline *s, double z){     //evaluates the spline
	int i=binsearchaspl(s->n,z,s);
	double h=z-s->x[i];
	return s->y[i]+h*(s->Ai[i]+h*(s->c[i]+h*s->d[i]));
}

double aspline_deriv(aspline *s, double z){     //evaluates d/dx of the spline
	int i=binsearchaspl(s->n,z,s);
	double h=z-s->x[i];
	return s->Ai[i]+2*h*s->c[i]+3*h*h*s->d[i];
}

double aspline_deriv2(aspline *s, double z){     //evaluates d/dx of the spline
	int i=binsearchaspl(s->n,z,s);
	double h=z-s->x[i];
	return 2*s->c[i]+6*h*s->d[i];
}

double aspline_integ(aspline* s, double z){//calculates the derivative of the aspline.
	int i = binsearchaspl(s->n,z,s);
	double sum=0;
	for(int k=0;k<i;k++){//integral for each interval is found and summed
		double xk=s->x[k], xk1=s->x[k+1], Ak=s->Ai[k], ck=s->c[k];
		double yk=s->y[k], dk=s->d[k], dx=xk1-xk;
		sum=sum+yk*dx+Ak*0.5*dx*dx + ck*dx*dx*dx/3+ dk*dx*dx*dx*dx/4;}
	double dxi=z-s->x[i];
	sum+= s->y[i]*dxi+s->Ai[i]*0.5*dxi*dxi + s->c[i]*dxi*dxi*dxi/3 + s->d[i]*dxi*dxi*dxi*dxi/4;
	return sum;
}

void aspline_free(aspline *s){ //free the allocated memory
	free(s->x); free(s->y); free(s->Ai); free(s->c); free(s->d); free(s);
}
