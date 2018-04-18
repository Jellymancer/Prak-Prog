#include<stdio.h>
#include<unistd.h> /* short comman-line options */
#include<getopt.h> /* long and short comman-line options */
#include<math.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_odeiv2.h>

int odefunc(double phi, const double y[], double f[], void *params){
double eps = *(double *) params;
f[0]=y[1];
f[1]=1-y[0]+eps*y[0]*y[0];
return GSL_SUCCESS;
}


int  main(){

/* assigment i */

double eps=0.01;
int dim=2;
    gsl_odeiv2_system sys = {odefunc, NULL, dim, &eps};
    gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

    int i;
    double phi0 = 0.0,  phif = 550.0; /* start and end of integration interval */
    double phi = phi0;
    double y[2] = {1.0,-0.5};  /* initial values */



    for (i = 1; i <= 1000; i++)
    {
        double phii = phi0 + i*(phif-phi0)/10000;
        int status = gsl_odeiv2_driver_apply (driver, &phi, phii, y);

        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;
        }
        printf("%.8e %.8e\n", phi, y[0]);
}
	gsl_odeiv2_driver_free(driver);
}
