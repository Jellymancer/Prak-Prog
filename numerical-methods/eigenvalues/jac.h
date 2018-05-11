
#ifndef HAVE_JAC_H
#define HAVE_JAC_H

#include<stdio.h>
#include<stdlib.h>
#define RND (double)rand()/RAND_MAX
#include<time.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
void printm(gsl_matrix*,FILE*);
void printv(gsl_vector*,FILE*);
int jac(gsl_matrix* , gsl_vector*,gsl_matrix*);
double* jac_onevalue(gsl_matrix* , gsl_vector*,gsl_matrix*,int,int);
double jac_eigvlbyeigvl(gsl_matrix* , gsl_vector*,gsl_matrix*);



#endif


