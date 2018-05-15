#include "integ.h"
#include <omp.h>

int main(){
//assignment 1: definite integrals
double sqr(double x,double c){return sqrt(x);}
double inversesqrt(double x,double c){return pow(sqrt(x),-1);}
double lnoversqrt(double x,double c){return log(x)/sqrt(x);}
double eq1(double x,double c){return 4*sqrt(1-pow((1-x),2));}

double gaussian(double x,double c){return exp(-x*x);}
double eq2(double x, double c){return 1/(x*x);};
double GSLgauss(double x, void* params){return exp(-x*x);}
//GSL functions are used for comparing with my integrator.
double GSLeq2(double x, void* params){return 1/(x*x);}


printf("Multiprocessing will be used to run the adaptive integrator. Two threads are made, one for assignment A"); 
printf(" and the other for assignment B.\n");

#pragma omp parallel sections
{
#pragma omp section
{
double x_start=0;
double x_end=1;
double abstol=0.00001;
double reltol=0.00001;
double result;
int *nrecs;
int n=0;
nrecs = &n;
result=integmain(sqr,x_start,x_end,abstol,reltol,nrecs);
printf("1:Assignment A: Recursive adaptive integration. Integrating the given four function and comparing the result to the listed value.\n\n");
printf("1:Numerical integration of sqrt(x) from 0 to 1 yields:  I=%f. It should be 2/3.\n",result);
printf("1:No. of evaluations: N=%d\n\n",*nrecs);

nrecs = &n;
result=integmain(inversesqrt,x_start,x_end,abstol,reltol,nrecs);
printf("1:Numerical integration of sqrt(x)⁻¹ from 0 to 1 yields:  I=%g. It should be 2.\n",result);
printf("1:No. of evaluations: N=%d\n\n",*nrecs);

nrecs = &n;
result=integmain(lnoversqrt,x_start,x_end,abstol,reltol,nrecs);
printf("1:Numerical integration of ln(x)/sqrt(x) from 0 to 1 yields: I= %g. It should be -4.\n",result);
printf("1:No. of evaluations: N=%d\n\n",*nrecs);

nrecs= &n;
result=integmain(eq1,x_start,x_end,abstol,reltol,nrecs);
printf("1:Numerical integration of 4*sqrt(1-(1-x)²) from 0 to 1 yields: I=%g. It should be %g.\n",result,M_PI);
printf("1:No. of evaluations: N=%d\n\n",*nrecs);
printf("Thread 1 closed.\n");
}

//Assignment 2  indefinite integrals.
#pragma omp section
{
double x_start, x_end;
double abstol=0.00001;
double reltol=0.00001;
double result;
int *nrecs;
int n=0;
nrecs= &n;

gsl_function f; f.function = GSLgauss;
int limit = 10000;
gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(limit);
// setting workspace for GSL integrator and performing the integration.
double GSLresult, GSLerr;
gsl_integration_qagi(&f,abstol,reltol,limit-1,workspace,&GSLresult,&GSLerr);
x_start=-INFINITY; x_end=INFINITY;

result=integmain(gaussian,x_start,x_end,abstol,reltol,nrecs);

printf("2:\n\nAssigment B: Infinite integrals.\nThe used routines are compared to the results of GSL routines.\n");
printf("2:Numerical integration of exp(-x²) from -Inf to Inf.\n I=%g. I_GSL= %g.\n",result,GSLresult);
printf("2:No. of evaluations: %d\n\n",*nrecs);

x_start=5; x_end=INFINITY;
gsl_function f2; f2.function = GSLeq2; GSLresult=0; GSLerr=0;
gsl_integration_qagiu(&f2,x_start,abstol,reltol,limit-1,workspace,&GSLresult,&GSLerr);
gsl_integration_workspace_free(workspace);
result=integmain(eq2,x_start,x_end,abstol,reltol,nrecs);
printf("2:Numerical integration of 1/x² from 5 to Inf yields.\n I=%g. I_GSL= %g.\n",result,GSLresult);
printf("2:No. of evaluations: N=%d\n\n",*nrecs);
printf("Thread 2 closed.\n");
}
}

return 0;
}
