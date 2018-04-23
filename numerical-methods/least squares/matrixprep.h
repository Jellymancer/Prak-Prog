#ifndef HAVE_MATRIXPREP_H
#define	HAVE_MATRIXPREP_H

#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<stdio.h>

double funs(int i, double x);
void matrix_generation (gsl_matrix* M,gsl_vector* b, double* x, double* y,double* dy);
void printm(gsl_matrix* M,FILE*);
void printv(gsl_vector* V,FILE*);

#endif
