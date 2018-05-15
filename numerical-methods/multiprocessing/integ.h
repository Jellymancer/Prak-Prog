

#ifndef HAVE_INTEG_H
#define HAVE_INTEG_H

#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_integration.h>
double integmain(double f(double,double), double, double , double, double,int*);



#endif




