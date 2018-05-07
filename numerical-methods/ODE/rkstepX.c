#include "functions.h"

void rkstepX(
	double x,                                  /* the current value of the variable */
	double step,                                  /* the step to be taken */
	gsl_vector* yx,                                /* the current value y(x) of the sought function */
	void f(double x, gsl_vector* yx, gsl_vector* dydx), /* the right-hand-side, dydt = f(x,y) */
	gsl_vector* yxh,                               /* output: y(x+h) */
	gsl_vector* interr)                                /* output: error estimate dx */

{
int n=yx->size; // size

gsl_vector* k0=gsl_vector_alloc(n); //paramters k_0, k_1/2 and k used for the Runge-Kutta midpoint method are allocated.
gsl_vector* yx1=gsl_vector_alloc(n);
gsl_vector* k12=gsl_vector_alloc(n);


//Now i find the value of k used for the step advancment and the integration error (err).
f(x,yx,k0);
gsl_vector_memcpy(interr,k0);

gsl_vector_scale(k0,step/2);
gsl_vector_memcpy(yx1,yx);
gsl_vector_add(yx1,k0);

f(x+step/2,yx1,k12);
gsl_vector_scale(k12,step);
gsl_vector_memcpy(yxh,yx);
gsl_vector_add(yxh,k12);

f(x+step/2,yx1,k12);
gsl_vector_sub(interr,k12);
gsl_vector_scale(interr,step/2);

gsl_vector_free(k0);
gsl_vector_free(yx1);
gsl_vector_free(k12);
}
