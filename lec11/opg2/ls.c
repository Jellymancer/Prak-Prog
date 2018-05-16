
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_multimin.h>
#include<assert.h>
#define MODEL(t) A*exp(-(t)/B)+C


/* least squares model and the experimental data. Data is given as parameters. */


struct exp_data {int n; double *t,*y,*e;};



double
my_f (const gsl_vector *v, void *params)
{


double  A = gsl_vector_get(v, 0);
double  B = gsl_vector_get(v, 1);
double  C = gsl_vector_get(v, 2);

  
  struct exp_data *p= (struct exp_data*) params;
  	int n=p->n;
	double *t=p->t;
	double *y=p->y;
	double *e=p->e;
	double sum=0;
	for(int i=0;i<n;i++) sum+=pow( (MODEL(t[i]) - y[i]) /e[i] ,2);
return sum;
}

int main (void)
{

double t[]= {0.47,1.41,2.36,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
double y[]= {5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
double e[]= {0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};
int n = sizeof(t)/sizeof(t[0]);

  size_t iter = 0;
  int status;
  int dim =3; /* dimension of the problem*/
  const gsl_multimin_fminimizer_type *T;
  gsl_multimin_fminimizer *s;
  

  struct exp_data params;
  params.n=n;
  params.t=t;
  params.y=y;
  params.e=e;

  gsl_vector *x;
  gsl_vector *xs;
  gsl_multimin_function my_func;
  my_func.n = dim;
  my_func.f = my_f;
  my_func.params=(void*)&params;

  /* Starting poionts and steppting values */
  x = gsl_vector_alloc (dim);
  gsl_vector_set (x, 0, 5.0);
  gsl_vector_set (x, 1, 7.0);
  gsl_vector_set (x, 2, 1.0);

  xs = gsl_vector_alloc(dim);
  gsl_vector_set(xs,0,2);
  gsl_vector_set(xs,1,2);
  gsl_vector_set(xs,2,2);



  T = gsl_multimin_fminimizer_nmsimplex2;
  s = gsl_multimin_fminimizer_alloc (T, dim);

  gsl_multimin_fminimizer_set (s, &my_func,x, xs);

      printf("iterations, A, B(lifetime), C,\n\n");

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);

      if (status)
        break;
      status = gsl_multimin_test_size (s->size, 1e-3);

      if (status == GSL_SUCCESS) printf ("Minimum found at:\n");

	printf ("%5d, %.5f, %.5f, %.5f\n", iter,
              gsl_vector_get (s->x, 0), 
              gsl_vector_get (s->x, 1), 
              gsl_vector_get (s->x, 2));

    }

  while (status == GSL_CONTINUE && iter < 1000);

  printf("\n\nThe optimization was performed on the function f(A,B,C)=A*exp(-t/B)+C\n");
  printf("Best lifetime (B) is estimated to be B=%.5f [unit]\n",gsl_vector_get (s->x, 1));
  printf("The fit is plotted in figure.\n");
plotdata(gsl_vector_get (s->x, 0),gsl_vector_get (s->x, 1),gsl_vector_get (s->x, 2));

  gsl_multimin_fminimizer_free (s);
  gsl_vector_free (x);
  gsl_vector_free (xs);
  return 0;
}




int plotdata(double A, double B, double C){

int i;
double j;
double t[]= {0.47,1.41,2.36,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
double y[]= {5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
double e[]= {0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};
int n = sizeof(t)/sizeof(t[0]);


FILE *fp;
fp = fopen("exp.dat","w");
for(i=0;i<=n;i++){
fprintf(fp, "%f %f %f\n", t[i],y[i],e[i]);
}
fclose(fp);

FILE *fp2;
fp2 = fopen("fit.dat","w");
for(j=0.0;j<=10;j+=0.1){
fprintf(fp2, "%f %f\n", j, A*exp(-(j)/B)+C);
}
fclose(fp2);



return 0;
}
