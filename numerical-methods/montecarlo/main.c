#include "mc.h"

int main(){
	int n=3, N=1000000; //allocating memory for variables. Set the number of random points (N) and the dimentsion of the integral (n).
	double val,err;
	gsl_vector* a=gsl_vector_alloc(n); gsl_vector* b=gsl_vector_alloc(n);
	double* res;
//Testing the algorythm.
	printf("\n\nAssignment 1: Some easy integrals are used for testing the algorythm.\n");
	printf("I first try to integrate a spherical volume element over a sphere using spherical coordinates.");
	printf("f=1, dV=r*r*sin(theta), r=[0;1] theta=[0;2*pi] phi=[0;pi].\n");
	gsl_vector_set_zero(a);
	gsl_vector_set(b,0,1); gsl_vector_set(b,1,M_PI);gsl_vector_set(b,2,2*M_PI);
	double f_sphere(gsl_vector* x){double r=gsl_vector_get(x,0); double theta=gsl_vector_get(x,1); return r*r*sin(theta);}
	res=plainmc(f_sphere,a,b,N);
	val=*res; err=*(res+1);
	printf("The result is I=%g. It should be 4/3*pi*r³=%g\nThe error is found to be: %g\n",val,4*M_PI/3,err);

	printf("\nNow integrating f=x*y*z over a cube in cartesian coordinates, dV=1, x=[0;1] y=[0;1] z=[0;1].\n");
	gsl_vector_set_zero(a);
	gsl_vector_set_all(b,1);
	double f_cube(gsl_vector* x){double xx=gsl_vector_get(x,0); double y=gsl_vector_get(x,1); double z=gsl_vector_get(x,2);
	return xx*y*z;}
	res=plainmc(f_cube,a,b,N);
	val=*res; err=*(res+1);
	printf("The result is I=%g. It should be 1/8=%g\nThe error is found to be: %g\n",val,0.125,err);


	printf("\nNow integrating the function given in assignment 1. f=[1-cos(x)cos(y)cos(z)]⁻¹ in (transformed) cartesian,");
	printf(" dV=1/pi, x=[0;pi] y=[0;pi] z=[0;pi]\n");
	gsl_vector_set_zero(a);
	gsl_vector_set_all(b,M_PI);
	double f_asg1(gsl_vector* x){double xx=gsl_vector_get(x,0); double y=gsl_vector_get(x,1); double z=gsl_vector_get(x,2);
	return pow(1-cos(xx)*cos(y)*cos(z),-1)*1/pow(M_PI,3);}
	res=plainmc(f_asg1,a,b,N);
	val=*res; err=*(res+1);
	printf("The result is I=%g. It should be I=%g\nThe error is found to be: %g\n",val,1.3932039296856768591842462603255,err);

//Checking error scaling
	printf("\nNow I check is the error of the MC integration method scales as 1/sqrt(N).\n");
	printf("This is done for the sphere integral by increasing N gradually. The result can be seen in figure 1.\n");

	FILE* errcheck; errcheck=fopen("errcheck.txt","w");
	for(int i=1;i<=1000;i++){
		gsl_vector_set_zero(a); gsl_vector_set_all(b,M_PI);
		double* res2=plainmc(f_cube,a,b,i); double err2 = *(res2+1);
		fprintf(errcheck,"%d %g\n",i,err2);
	}
	fclose(errcheck);

return 0;
}
