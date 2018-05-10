#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include "neural_network.h"

double activation_function(double x)
{
	return x*exp(-x*x);
}

double fit_func(double x)
{
	return cos(5*x-1)*exp(-x*x);
}

int main()
{
	int n = 20;
	neural_network* NN = neural_network_alloc(n, activation_function);
	double a = -1, b = 1;
	int nx = 10;
	gsl_vector* x_list = gsl_vector_alloc(nx);
	gsl_vector* y_list = gsl_vector_alloc(nx);

	for (int i = 0; i < nx; i++)
	{
		double x = a +  (b-a)*i/(nx-1);
		double f = fit_func(x);
		gsl_vector_set(x_list, i, x);
		gsl_vector_set(y_list, i, f);
	}

	for (int i = 0; i < NN->size; i++)
	{
		gsl_vector_set(NN->data, 0*NN->size + i, a + (b-a)*i/(NN->size - 1));
		gsl_vector_set(NN->data, 1*NN->size + i, 1);
		gsl_vector_set(NN->data, 2*NN->size + i, 1);
	}

	neural_network_trainer(NN, x_list, y_list);

	for (int i = 0; i < x_list->size ; i++)
	{
		double x = gsl_vector_get(x_list, i);
		double f = gsl_vector_get(y_list, i);
		printf("%g %g \n", x, f);
	}
	printf("\n\n");

	double dz = 1.0/(2*64);
	
	for(double z = a; z <= b; z += dz){
	
		double y = neural_network_feed_forward(NN,z);

		printf("%g %g\n", z, y);
	}

	neural_network_free(NN);
	gsl_vector_free(x_list);
	gsl_vector_free(y_list);

	return 0;
}

