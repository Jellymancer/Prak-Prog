#ifndef HAVE_NEURONS_H
#define HAVE_NEURONS_H
#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_multimin.h>
#include<assert.h>
typedef struct {
	int n;
	double(*fun)(double);
	gsl_vector* data;
	} neurons;
neurons* neurons_alloc(int n,double(*fun)(double));

void neurons_free(neurons* nw);
double neurons_feed_forward(neurons* nw, double x);
void neurons_train(neurons* nw,gsl_vector* xlist, gsl_vector* ylist);
#endif

