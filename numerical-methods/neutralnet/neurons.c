#include "neurons.h"

neurons* neurons_alloc(int n, double(*fun)(double)){
	neurons* nw =malloc(sizeof(neurons));
	nw->n=n;
	nw->fun=fun;
	nw->data=gsl_vector_alloc(3*n);
	return nw;
}

void neurons_free(neurons* nw){
	gsl_vector_free(nw->data);
	free(nw);}

double neurons_feed_forward(neurons* nw, double x){
	double s=0;
	for(int i=0;i<nw->n;i++){
		double a=gsl_vector_get(nw->data,0*nw->n+i); /*vector has a,b,w */
		double b=gsl_vector_get(nw->data,1*nw->n+i);
		assert(b != 0);
		double w=gsl_vector_get(nw->data,2*nw->n+i);
		s+=nw->fun((x+a)/b)*w; //input signal is transformed by every neuron
	}
	assert(isnan(s)==0 && isinf(s)==0);
	return s;
}

void neurons_train(neurons* nw,gsl_vector* vx, gsl_vector* vy){
//This function trains the network by comparing the network-output values to values of the tabulated function.
//GSL library minimalization routine is then employed for optimizing the weigths.
	double delta (gsl_vector* p, void* params) {
	gsl_vector_memcpy(nw->data,p);
	//p are the parameters to be optimized by the network.
	double s=0;
	for(int i=0;i<vx->size;i++)
		{
		double x=gsl_vector_get(vx,i); // input x, output f and parameters y;
		double fy=gsl_vector_get(vy,i);
		double y=neurons_feed_forward(nw,x);
		s+=fabs(y-fy); // calculate the difference between the calculated and correct result.
		}
		assert(isnan(s)==0 && isinf(s)==0);
		return s/(vx->size);
	}

	gsl_vector* p=gsl_vector_alloc(nw->data->size); //allocating memory and setting start values for the GSL routine used for training the network.
	gsl_vector_memcpy(p,nw->data);
	gsl_vector* stp=gsl_vector_alloc(nw->data->size); //stepping values.
	gsl_vector_set_all(stp,0.01);

	gsl_multimin_function fGSL;
	fGSL.f = delta;
	fGSL.n = nw->data->size;
	fGSL.params = NULL;
	const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s =gsl_multimin_fminimizer_alloc(T,nw->data->size);
	gsl_multimin_fminimizer_set (s, &fGSL,p,stp);
  int iter=0, status1, status2;
  do
    {
      iter++;
      status1 = gsl_multimin_fminimizer_iterate(s);

      if (status1 !=0){printf("GSL MINIMIZER ERROR\n"); break;}
      status2 = gsl_multimin_test_size (s->size, 1e-6);
    }
  while (status2 == GSL_CONTINUE && iter < 10000);

	gsl_vector_memcpy(nw->data,s->x);
	gsl_vector_free(p);
	gsl_vector_free(stp);
	gsl_multimin_fminimizer_free(s);
}

