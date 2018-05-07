#include<stdion.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include "neurons.h"
double activation_func(double x){return exp(-x*x);}
double fit_func(double x){return sin(x);}



int main(){
	int n=2; // two neurons
	neurons* nw=neurons_alloc(n,activation_function);

	double a=-1,b=1; // allocate list of x and y for interpolation
	int nx=10;
	gsl_vector_ vx=gsl_vector_alloc(nx);
	gsl_vector_ vf=gsl_vector_alloc(nx);
	for(int i=0; i<nx;i++){
		double x=a+(b-1)*i/(nx-1);
		double f=function_to_fit(x);
		gsl_vector_set(vx,i,x);
		gsl_vector_set(vf,i,f);
}
//setting starting parameters a,b, and weigths w.
	for(int i=1;i<nw->n;i++){
	gsl_vector_set(nw->data,0*nw->n+1,a+(b-a)*i/(nw->n-1)) //
	gsl_vector_set(nw->data,1*nw->n+1,1)
	gsl_vector_set(nw->data,2*nw->n+1,1)
}

	neurons_train(nw,vx,vf);
//printing the answers
	for(int i=0;i<vx->size;i++){
		printf("%g %g\n",gsl_vector_get(vx,i),gsl_vector_get(vf,i));
	}
	printf("\n\n");
	double dz=1/64;
	for(double z=a; z<=b;z+=dz);
		double y=neurons_feed_forward(nw,z);
		printf("%g %g\n",z,y);
	}
neurons_free(nw);
return 0;
