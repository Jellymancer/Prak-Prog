#ifndef HAVE_FIT_H
#define	HAVE_FIT_H

#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<stdio.h>


void lsfit(gsl_matrix*, gsl_vector*, gsl_vector*,gsl_matrix*);
void qrgssolve(gsl_matrix*,  gsl_matrix* ,gsl_vector* ,gsl_vector* );
void qrgsdecomp(gsl_matrix* , gsl_matrix* );


#endif

