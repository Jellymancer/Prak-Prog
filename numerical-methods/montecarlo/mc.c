#include "mc.h"
#define RND ((double) rand ()/RAND_MAX)

void randabs(gsl_vector* a, gsl_vector* b, gsl_vector* x){
	for(int i=0;i<x->size;i++){
		double ai = gsl_vector_get(a,i), bi = gsl_vector_get(b,i);
		gsl_vector_set(x,i,ai+RND*(bi-ai));
	}
}


double* plainmc(double f(gsl_vector* x),gsl_vector* a,gsl_vector* b,int N){
	int dim=a->size;
	gsl_vector* x=gsl_vector_alloc(dim);
	double V=1; for(int i=0;i<dim;i++) V*=gsl_vector_get(b,i)-gsl_vector_get(a,i);

	double sum=0,sum2=0;
	for(int i=0;i<N;i++){randabs(a,b,x); double fx=f(x); sum +=fx; sum2 += pow(fx,2);}

	double mean = sum/N;
	double var = sum2/N - mean*mean;
	static double result[2];
	result[0]=mean*V; result[1]=sqrt(var/N)*V;
	gsl_vector_free(x);
	return result;
}



