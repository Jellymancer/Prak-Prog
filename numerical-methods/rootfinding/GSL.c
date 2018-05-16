#include "func.h"

int GSLrootf(int GSLf(const gsl_vector *v, void *params, gsl_vector *f),gsl_vector* v,double eps){

size_t iter = 0;
int status;
int dim = v->size;
const gsl_multiroot_fsolver_type *T;
gsl_multiroot_fsolver *s;

gsl_multiroot_function my_func;
my_func.n = dim;
my_func.f = GSLf;
my_func.params=NULL;


T = gsl_multiroot_fsolver_hybrids;
s = gsl_multiroot_fsolver_alloc(T,dim);

gsl_multiroot_fsolver_set(s,&my_func,v);

do
    {
      iter++;
      status=gsl_multiroot_fsolver_iterate (s);
	if (status) break;
      status = gsl_multiroot_test_residual(s->f, eps);
    }
  while (status == GSL_CONTINUE && iter < 10000);

//printv(s->x,stdout); //For checking end value.
//printf("%zu\n\n\n",iter);
gsl_multiroot_fsolver_free(s);

return iter;
}
