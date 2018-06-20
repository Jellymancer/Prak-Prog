#ifndef HAVE_INTERPOL_H
#define HAVE_INTERPOL_H


#include<stdio.h>
#include<stdlib.h>
#define RND (double)rand()/RAND_MAX
#include<math.h>
#include<assert.h>


typedef struct {int n; double *x,*y,*b,*c,*d;} cspline;
cspline *cspline_alloc(int n, double *x, double *y);
double cspline_eval(cspline *s, double z);
double cspline_deriv(cspline *s, double z);
double cspline_deriv2(cspline *s, double z);
double cspline_integ(cspline *s, double z);

typedef struct {int n; double *x,*y,*b,*c,*d;} qspline;
qspline *qspline_alloc(int n, double *x, double *y);
double qspline_eval(qspline *s, double z);
double qspline_deriv(qspline *s, double z);
double qspline_deriv2(qspline *s, double z);
double qspline_integ(qspline *s, double z);

typedef struct {int n; double *x,*y,*df,*b,*c,*d;} ospline;
ospline *ospline_alloc(int n, double *x, double *y, double *df);
double ospline_eval(ospline *s, double z);
double ospline_deriv(ospline *s, double z);
double ospline_deriv2(ospline *s, double z);
double ospline_integ(ospline *s, double z);

typedef struct {int n; double *x, *y, *Ai, *c, *d;} aspline;
aspline* aspline_alloc(int n,double* x,double* y);
double aspline_eval(aspline *s, double z);
double aspline_deriv(aspline *s, double z);
double aspline_deriv2(aspline *s, double z);
double aspline_integ(aspline *s, double z);
void aspline_free(aspline *s);

int binsearchaspl(int n, double z, aspline* s);
int binsearchcspl(int n, double z, cspline* s);
int binsearchqspl(int n, double z, qspline* s);
int binsearchospl(int n, double z, ospline* s);
#endif
