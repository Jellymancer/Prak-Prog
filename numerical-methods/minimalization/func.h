
#ifndef HAVE_FUNC_H
#define HAVE_FUNC_H

#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
void qrgsdecomp(gsl_matrix *E,gsl_matrix *W);

void qrgssolve(gsl_matrix* Q,  gsl_matrix* R,gsl_vector* b,gsl_vector* c);

double fRos(gsl_vector* p,gsl_vector* fx);

void Jacobi_fRos(gsl_vector* p, gsl_matrix* J);

double fHim(gsl_vector* p, gsl_vector* fx);

void Jacobi_fHim(gsl_vector* p, gsl_matrix* J);

double fRad(gsl_vector* p,gsl_vector* fx);

int newton(double f(gsl_vector* x,gsl_vector* fx), void Hes(gsl_vector* p, gsl_matrix* H),gsl_vector*, double, double);

int newton_broy(double f(gsl_vector* x,gsl_vector* dfx),gsl_vector*, double);

void printm(gsl_matrix* M,FILE*);
void printv(gsl_vector* V,FILE*);


#endif
