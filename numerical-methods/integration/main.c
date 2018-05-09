#include "integ.h"
int main(){
//assignment 1: definite integrals
double sqr(double x,double c){return sqrt(x);}
double inversesqrt(double x,double c){return pow(sqrt(x),-1);}
double lnoversqrt(double x,double c){return log(x)/sqrt(x);}
double eq1(double x,double c){return 4*sqrt(1-pow((1-x),2));}

double x_start=0;
double x_end=1;
double abstol=0.0001;
double reltol=0.0001;
double result;
int *nrecs;
int n=0;

nrecs = &n;
result=integmain(sqr,x_start,x_end,abstol,reltol,nrecs);
printf("Assignment 1: Recursive adaptive integration. Integrating the following four functions.\n");
printf("Numerical integration of sqrt(x) from 0 to 1 yields:  I=%f. It should be 2/3.\n",result);
printf("No. of evaluations: N=%d\n\n",*nrecs);

nrecs = &n;
result=integmain(inversesqrt,x_start,x_end,abstol,reltol,nrecs);
printf("Numerical integration of sqrt(x)⁻¹ from 0 to 1 yields:  I=%g. It should be 2.\n",result);
printf("No. of evaluations: N=%d\n\n",*nrecs);

nrecs = &n;
result=integmain(lnoversqrt,x_start,x_end,abstol,reltol,nrecs);
printf("Numerical integration of ln(x)/sqrt(x) from 0 to 1 yields: I= %g. It should be -4.\n",result);
printf("No. of evaluations: N=%d\n\n",*nrecs);

nrecs= &n;
result=integmain(eq1,x_start,x_end,abstol,reltol,nrecs);
printf("Numerical integration of 4*sqrt(1-(1-x)²) from 0 to 1 yields: I=%g. It should be %g.\n",result,M_PI);
printf("No. of evaluations: N=%d\n\n",*nrecs);
printf("\n\nAssigment 2: Infinite integrals.\nThe used routines are compared to the results of GSL routines.\n");


//Assignment 2  indefinite integrals.
nrecs= &n;
double gaussian(double x,double c){return exp(-x*x);}
double eq2(double x, double c){return 1/(x*x);};

double GSLgauss(double x, void* params){return exp(-x*x);}
//GSL integration of the gaussian used to compare.
double GSLeq2(double x, void* params){return 1/(x*x);}

gsl_function f; f.function = GSLgauss;
int limit = 10000;
gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(limit);
double GSLresult, GSLerr;
gsl_integration_qagi(&f,abstol,reltol,limit-1,workspace,&GSLresult,&GSLerr);
x_start=-INFINITY; x_end=INFINITY;
result=integmain(gaussian,x_start,x_end,abstol,reltol,nrecs);
printf("Numerical integration of exp(-x²) from -Inf to Inf.\n I=%g. I_GSL= %g.\n",result,GSLresult);
printf("No. of evaluations: %d\n\n",*nrecs);

x_start=5; x_end=INFINITY;
gsl_function f2; f2.function = GSLeq2; GSLresult=0; GSLerr=0;
gsl_integration_qagiu(&f2,x_start,abstol,reltol,limit-1,workspace,&GSLresult,&GSLerr);
gsl_integration_workspace_free(workspace);
result=integmain(eq2,x_start,x_end,abstol,reltol,nrecs);
printf("Numerical integration of 1/x² from 5 to Inf yields.\n I=%g. I_GSL= %g.\n",result,GSLresult);
printf("No. of evaluations: N=%d\n\n",*nrecs);
return 0;
}
