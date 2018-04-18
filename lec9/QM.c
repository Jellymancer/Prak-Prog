#include<math.h>
#include<stdio.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>


double norm(double x, void* params){
double z=*(double*)params;
int p=2;
return exp((-1*z*pow(x,p)));
}


double hamil(double x,void* params){
double z=*(double*)params;
int p=2;
return ( (-1*pow(z,p)*pow(x,p))/2  + z/2 + pow(x,p)/2)*exp(-1*z*pow(x,p));
}



double integration(double z){
	gsl_function f;
	f.function = norm;
	f.params = (void*)&z;

	int limit = 100;
	gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(limit);
	double epsabs=0.000001, epsrel=0.000001, resultnorm, errnorm, resulthamil, errhamil;
	int status1=gsl_integration_qagi(&f,epsabs,epsrel,limit-1,workspace,&resultnorm,&errnorm);
	gsl_integration_workspace_free(workspace);

	gsl_function g;
        g.function = hamil;
        g.params = (void*)&z;
        gsl_integration_workspace* workspaceb = gsl_integration_workspace_alloc(limit);
        int status2=gsl_integration_qagi(&g,epsabs,epsrel,limit-1,workspaceb,&resulthamil,&errhamil);
        gsl_integration_workspace_free(workspaceb);

	if(status1!=GSL_SUCCESS & status2!=GSL_SUCCESS) return NAN;
	else return resulthamil/resultnorm;
}


int main(){
double a=0.01, b=10, dx=0.01;
double x;

for(x=a;x<=b;x+=dx)
	printf("%g %g\n",x,integration(x));


return 0;
}

