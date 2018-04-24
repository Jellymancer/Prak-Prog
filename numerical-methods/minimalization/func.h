
#ifndef HAVE_FUNC_H
#define HAVE_FUNC_H

#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<stdio.h>
#include<math.h>
void fRos(gsl_vector* p,gsl_vector* fx);

void Jacobi_fRos(gsl_vector* p, gsl_matrix* J);

void fSys(gsl_vector* p, gsl_vector* fx);

void Jacobi_fSys(gsl_vector* p, gsl_matrix* J);

void fHim(gsl_vector* p, gsl_vector* fx);

void Jacobi_fHim(gsl_vector* p, gsl_matrix* J);

int newton(void f(gsl_vector* x,gsl_vector* fx), void Jacobi(gsl_vector* p, gsl_matrix* J),gsl_vector* x, double dx, double eps, int numericJac);

#endif
