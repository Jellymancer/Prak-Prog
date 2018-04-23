#ifndef HAVE_MATRIXPREP_H
#define HAVE_MATRIXPREP_H

#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<stdio.h>

void printm(gsl_matrix* M,FILE*);
void printv(gsl_vector* V,FILE*);

#endif






