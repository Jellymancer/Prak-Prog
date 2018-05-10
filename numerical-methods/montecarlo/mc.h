
#ifndef HAVE_MC_H
#define HAVE_MC_H

#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<stdlib.h>

double* plainmc(double f(gsl_vector* x),gsl_vector*,gsl_vector*, int);



#endif


