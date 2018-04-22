
#ifndef HAVE_INTERPOL_H
#define HAVE_INTERPOL_H


#include<stdio.h>
#include<stdlib.h>
#define RND (double)rand()/RAND_MAX
#include<math.h>


double linterp(int,double*,double*,double);
double linterpint(int,double*,double*,double);
typedef struct {int n; double *x, *y, *b, *c;} qspline;
qspline* qspline_alloc(int n,double* x,double* y);
double qspline_eval(qspline *s, double z);
void qspline_free(qspline *s);
double qspline_deriv(qspline * s, double z);
double qspline_intergral(qspline *q, double z);

#endif
