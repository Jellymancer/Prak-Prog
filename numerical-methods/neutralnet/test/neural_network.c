#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>

typedef struct {int size; double (*function)(double); gsl_vector* data;} neural_network;

neural_network* neural_network_alloc(int number_of_hidden_neurons, double(*activation_function)(double)) 
{
	neural_network* NN = malloc(sizeof(neural_network*));
	NN->size = number_of_hidden_neurons;
	NN->function = activation_function;
	NN->data = gsl_vector_alloc(3*number_of_hidden_neurons);

	return NN;		
}

void neural_network_free(neural_network* NN) 
{
	gsl_vector_free(NN->data);
	free(NN);
}

double neural_network_feed_forward(neural_network* NN, double x)
{
	double s = 0;
	for (int i = 0; i < NN->size; i++) {

		double a = gsl_vector_get(NN->data, 0*NN->size + i);
		double b = gsl_vector_get(NN->data, 1*NN->size + i);
		double w = gsl_vector_get(NN->data, 2*NN->size + i);
		s += NN->function((x+a)/b)*w;
	}

	return s;
}

void neural_network_trainer(neural_network* NN, gsl_vector* xlist, gsl_vector* ylist)
{

	double delta(gsl_vector* p)
	{
		gsl_vector_memcpy(NN->data, p);
		double s = 0;
		for (int i = 0; i < xlist->size; i++)
		{
			double x = gsl_vector_get(xlist, i);
			double f = gsl_vector_get(ylist, i);
			double y = neural_network_feed_forward(NN, x);
			s += fabs(y - f);
		}
		return s/xlist->size;
	}

	gsl_vector* p = gsl_vector_alloc(NN->data->size);
	gsl_vector* step_size = gsl_vector_alloc(NN->data->size);
	gsl_vector_memcpy(p, NN->data); 
	gsl_vector_set_all(step_size, 0.1);

	gsl_multimin_function F;
	F.f = delta;
	F.n = p->size;
	F.params = NULL;

	gsl_multimin_fminimizer *s = 
	gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, F.n);		
	gsl_multimin_fminimizer_set (s, &F, p, step_size);
	
	int iter = 0, status;
	do{
		iter++;
		int iteration_status = gsl_multimin_fminimizer_iterate(s);
		if(iteration_status != 0)
		{
			fprintf(stderr, "unable to improve\n");
			break;
		}
		double acc = 0.001;
		status = gsl_multimin_test_size(s->size, acc);
		if(status == GSL_SUCCESS) 
		{	
			fprintf(stderr, "neural network converged in %i iterations\n", iter);
		}

	}while( status == GSL_CONTINUE && iter < 1e+6);

	gsl_vector_memcpy(NN->data, s->x); 

	
	fprintf(stderr, "neural network parameters\n");
	fprintf(stderr, "a \t b \t w\n");
	for (int i = 0; i < NN->size; i++)
	{
		fprintf(stderr, "%g \t %g \t %g \n", 
			gsl_vector_get(NN->data, 0*NN->size + i),
			gsl_vector_get(NN->data, 1*NN->size + i),
			gsl_vector_get(NN->data, 2*NN->size + i) );
	}

	
	gsl_vector_free(p);
	gsl_vector_free(step_size);	
	gsl_multimin_fminimizer_free(s);
}

