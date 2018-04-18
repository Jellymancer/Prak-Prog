
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_matrix.h>

int func(double x, const double y[], double f[], void *params){
f[0]=y[0]*(1-y[0]);
return GSL_SUCCESS;
}
/* function above determines the diff equation to be integrated, x is the parameter, y is the function and f is the derivative dydx) */

int main(void)
{
int dim=1;
gsl_odeiv2_system sys={func, NULL ,dim, NULL};
gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,1e-6,1e-6,0.0);
/* inicialization of ODE function in GSL, NULL values here as this is only a first order problem */

int i;
double x0=0.0, xf=5.0; /* start, end and stepsize of integration interval */
double x=x0;
double y[1]={0.5}; /* initial value*/
for(i=1;i<=100;i++){

double xi=x0+i*(xf-x0)/100;
int status = gsl_odeiv2_driver_apply(d, &x,xi,y);

if(status!=GSL_SUCCESS)
{printf("ERROR, return value=%d",status);
break;}


printf("%.8e %.8e\n",x,y[0]);
}

gsl_odeiv2_driver_free(d);
return 0;
}

