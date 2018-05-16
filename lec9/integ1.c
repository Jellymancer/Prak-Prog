#include<math.h>
#include<stdio.h>
#include<gsl/gsl_integration.h>



double integrand(double x, void* params){
double z=*(double*)params;
return log(x)/sqrt(x);
}


double funky(double z){
	gsl_function f;
	f.function = integrand;
	f.params = (void*)&z;

	int limit = 100;
	gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(limit);
	double a=0, b=1, acc=100, eps=100, result, err;
	gsl_integration_qags(&f,a,b,acc,eps,limit-1,workspace,&result,&err);
	gsl_integration_workspace_free(workspace);

	return result;
}


int main(){
int x=1;
printf("\n\nAssignment 1: Performing the integration yields the value %g\n",funky(x));
printf("Assignment 2 is plotted in the attched figure\n");
return 0;
}

