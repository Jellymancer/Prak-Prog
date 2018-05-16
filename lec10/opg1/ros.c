#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>
#include<stdlib.h>
#include<stdio.h>

struct rparams
{
double a;
double b;
};
int equation(const gsl_vector* v, void* params, gsl_vector* f){
	double a = ((struct rparams*) params)->a;
	double b = ((struct rparams*) params)->b;
        const double v0=gsl_vector_get(v,0);
        const double v1=gsl_vector_get(v,1);
        const double f0= a*(v1-v0);
	const double f1= b*v0*v0*v0-2*v0*(a*v1-1)-2;
	gsl_vector_set(f,0,f0);
        gsl_vector_set(f,1,f1);
return GSL_SUCCESS;
}

int main(void){
	int status;
        const gsl_multiroot_fsolver_type* T =gsl_multiroot_fsolver_hybrids;
        gsl_multiroot_fsolver* S = gsl_multiroot_fsolver_alloc(T,2);
        gsl_multiroot_function F;
	int a=1;
	int b=100;
	struct rparams p= {a , b};
        F.f=&equation; /* Struct F containing the function, its dimension and parameters */
        F.n=2;
        F.params=&p;
        gsl_vector *x = gsl_vector_alloc(2);
        gsl_vector_set(x,0,20);
        gsl_vector_set(x,1,20);
        gsl_multiroot_fsolver_set(S,&F,x);
	int iter=0; /*number of iterations*/

	print_state (iter,S);


	    do {
      	iter++;
      	status = gsl_multiroot_fsolver_iterate (S);

      print_state (iter, S);

      if (status)   /* check if solver is stuck */
        break;

      status =
        gsl_multiroot_test_residual (S->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s\n", gsl_strerror (status));
  printf ("The extremum of the Rosenbrock function is found at x=%.3f y=%.3f\n",gsl_vector_get(S->x,0),gsl_vector_get(S->x,1));
  printf ("It took %3u iterations to find the extremum\n", iter);
  gsl_multiroot_fsolver_free (S);
  gsl_vector_free (x);
  return 0;
}



int
print_state (size_t iter, gsl_multiroot_fsolver* S)
{
  printf ("iterations = %3u x = %.3f y=%.3f "
          "dfdy = %.3e dfdx=%.3e\n",
          iter,
          gsl_vector_get (S->x, 0),
          gsl_vector_get (S->x, 1),
          gsl_vector_get (S->f, 0),
          gsl_vector_get (S->f, 1));
}

