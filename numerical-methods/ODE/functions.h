
#ifndef HAVE_FUNCTIONS_H
#define HAVE_FUNCTIONS_H

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_sf_airy.h>
void rkstepX(double , double , gsl_vector*, void f(double, gsl_vector*, gsl_vector*) ,gsl_vector* , gsl_vector* );

int driver(double ,double ,double ,gsl_vector*,double,double,int,
void rkstepX(double , double , gsl_vector*,
void f(double ,gsl_vector* ,gsl_vector* ),
gsl_vector* , gsl_vector*),
void f(double ,gsl_vector *,gsl_vector *));

int driverwpath(double ,double ,double ,gsl_vector*,double,double,int,
void rkstepX(double , double , gsl_vector*,
void f(double ,gsl_vector* ,gsl_vector* ),
gsl_vector* , gsl_vector*),
void f(double ,gsl_vector *,gsl_vector *),
gsl_matrix*, int);




#endif


