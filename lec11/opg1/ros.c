#include<stdio.h>
#include<math.h>
#include<gsl/gsl_multimin.h>
#include<assert.h>


/* Rosenbrock function with two parameters a=1=p[0] and b=100=p[1] */

double
my_f (const gsl_vector *v, void *params)
{
  double x, y;
  double *p = (double *)params;
  
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);
 
  return pow((p[0]-x),2) + p[1] * pow(y-pow(x,2),2);
}

/* The gradient of f, df = (df/dx, df/dy). */
void 
my_df (const gsl_vector *v, void *params, 
       gsl_vector *df)
{
  double x, y;
  double *p = (double *)params;
  
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);
 
  gsl_vector_set(df, 0, 4*p[1]*x*x*x - 2*x*(2*p[1]*y-1) - 2*p[0]);
  gsl_vector_set(df, 1, 2*p[1]*y - 2*p[1]*x*x);
}

/* Compute both f and df together. */
void 
my_fdf (const gsl_vector *x, void *params, 
        double *f, gsl_vector *df) 
{
  *f = my_f(x, params); 
  my_df(x, params, df);
}

int main (void)
{
  size_t iter = 0;
  int status;
  int dim =2; /* dimension of the problem*/
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  /* Parameter definition, here a=1 and b=100 */
  double par[5] = { 1.0, 100.0 };

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.n = dim;
  my_func.f = my_f;
  my_func.df = my_df;
  my_func.fdf = my_fdf;
  my_func.params = par;

  /* Starting point, x = (5,7) */
  x = gsl_vector_alloc (2);
  gsl_vector_set (x, 0, 5.0);
  gsl_vector_set (x, 1, 7.0);

  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, 2);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);

      printf("iterations, x, y, f(x,y)\n\n");

  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      if (status)
        break;

      status = gsl_multimin_test_gradient (s->gradient, 1e-3);

      if (status == GSL_SUCCESS)
        printf ("Minimum found at:\n");


      printf ("%5d, %.5f, %.5f, %10.5f\n", iter,
              gsl_vector_get (s->x, 0), 
              gsl_vector_get (s->x, 1), 
              s->f);

    }
  while (status == GSL_CONTINUE && iter < 1000);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);

  return 0;
}

