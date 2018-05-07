/*typedef struct {
        int n;
        doble(*f)(double);
        gsl_vector* data;
        } neurons;
neurons* neurons_alloc(int n,double(*f)(double)):

void neurons_free(neurons* nw);
double neurons_feed_forward(neurons* nw, double x);
void neurons_train(neurons* nw,gsl_vector* xlist, gsl_vector* ylist);
*/
#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<math.h>
#include "neurons.h"

neurons* neurons_alloc(int n, double(*f)(double)){
	neurons* nw =malloc(sizeof(neurons));
	nw->n=n;
	nw->f=f;
	nw->gsl_vector_alloc(3*n);
	return nw;
}

void neurons_free(neurons* nw){
	gsl_vector_free(nw->data;
	free(nw);}

double neurons_feed_forward(neurons* nw, double x){
	double s=0;
	for(int i=0;i<nw->n;i++){
		double a=gsl_vector_get(nw->data,0*nw->n+i,) /*vector has a,b,w */
		double b=gsl_vector_get(nw->data,1*nw->n+i)
		double w=gsl_vector_get(nw->data,2*nw->n+i)
		s+=nw->f((x+a)/b);
	}
// needs fixing
	retun s;

int qnewton(double deviation(gsl_vector* x), gsl_vector* x, double eps);

void neurons_train(neurons* nw,gsl_vector* vx, gsl_vector* vy){
// For this use GSL amoeba
//This function trains the network.
	double delta(gsl_vector* p){
//p are the parameters to ba adjusted by the network.
}	double s=0;
	for(int i=1;i<vx->size;i++){
		double x=gsl_vector_get(vx,i); // input x, output f and parameters y;
		double f=gsl_Vector_getvf,i);
		double y=neurons_feed_forward(nw,x);
		s+=fabs(y-f); // calculate the difference between the calculated and correct result.
	}
	return s/vx->size

	gsl_vector* p=gsl_vector_alloc(nw->data->size);
	gsl_vector_memcpy(p,nw->data);
	qnewton(...) //this is to be replaced with a different routine gsl prolly.
	gsl_vector_memcpy(nw->data,p);
	gsl_vector_free(p);
}

