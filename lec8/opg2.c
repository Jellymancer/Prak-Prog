#include<stdio.h>
#include<unistd.h> /* short comman-line options */
#include<getopt.h> /* long and short comman-line options */
#include<math.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_odeiv2.h>

int odefunc(double x, const double y[], double f[], void *params){
double eps = *(double *) params;
f[0]=y[1];
f[1]=1-y[0]+eps*y[0]*y[0];
return GSL_SUCCESS;
}


int  main(){

/* assigment i */

double eps=0.0;
int dim=2;
    gsl_odeiv2_system sys = {odefunc, NULL, dim, &eps};
    gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

    int i;
    double x0 = 0.0,  xf = 80.0; /* start and end of integration interval */
    double x = x0;
    double y[2] = {1.0, 0.0};  /* initial values */



    for (i = 1; i <= 100; i++)
    {
        double xi = x0 + i * (xf-x0) / 1000.0;
        int status = gsl_odeiv2_driver_apply (driver, &x, xi, y);

        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;
        }
        printf("%f %f\n", x, y[0]);
  }
	gsl_odeiv2_driver_free(driver);

}

