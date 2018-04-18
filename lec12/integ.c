#include<stdio.h>
#include<unistd.h> /* short comman-line options */
#include<getopt.h> /* long and short comman-line options */
#include<math.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_odeiv2.h>

int odefunc(double x, const double y[], double f[], void *params){
f[0]=(2/(sqrt(M_PI))) * exp(-1*pow(x,2));
return GSL_SUCCESS;
}


int  main(){

double x0;
double xf;
double dx;

printf( "The solving algorythm reuqires three inputs. a = start of integration, b = end of integration, dx = integration step. ");
printf("\n The integration should go from negative value to a positive one.");
printf("\n Please input the parameters (format: a b dx): ");
scanf("%lf %lf %lf", &x0, &xf, &dx);
printf("You entered: a=%lf, b=%lf, dx=%lf\n", x0, xf, dx);

int dim=1;
    gsl_odeiv2_system sys = {odefunc, NULL, dim, NULL};
    gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

    double i;
  /*  double x0 = -8.0,  xf = 8.0;  start and end of integration interval */
    double x = x0;
    double y[1] = {-1.0};  /* initial values */


	FILE *fg;
        fg = fopen("data.txt","w");

    for (i = x0; i<=xf+dx; i+=dx)
    {
        int status = gsl_odeiv2_driver_apply (driver, &x, i, y);

        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;
        }

        fprintf(fg,"%f %f\n", x, y[0]);
}
	fclose(fg); 

	gsl_odeiv2_driver_free(driver);


return 0;
}

