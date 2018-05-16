
#ifndef HAVE_FUNC_H
#define HAVE_FUNC_H

#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_multiroots.h>

void fRos(gsl_vector* p,gsl_vector* fx);

void Jacobi_fRos(gsl_vector* p, gsl_matrix* J);

void fSys(gsl_vector* p, gsl_vector* fx);

void Jacobi_fSys(gsl_vector* p, gsl_matrix* J);

void fHim(gsl_vector* p, gsl_vector* fx);

void Jacobi_fHim(gsl_vector* p, gsl_matrix* J);

int GSLfRos(const gsl_vector*, void *params,gsl_vector*);

int GSLfSys(const gsl_vector*, void *params,gsl_vector*);

int GSLfHim(const gsl_vector* p, void *params,gsl_vector* f);

int newton(void f(gsl_vector* x,gsl_vector* fx), void Jacobi(gsl_vector* p, gsl_matrix* J),gsl_vector* x, double dx, double eps, int numericJac);

int GSLrootf(int GSLf(const gsl_vector *v, void *params, gsl_vector *f),gsl_vector* startv,double eps);

#endif
