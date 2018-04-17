
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>
#include<stdio.h>
#include<math.h>


int root_equation(const gsl_vector* v, void* params, gsl_vector* f){
	double x=gsl_vector_get(v,0);
	double z=*(double*)params;
	gsl_vector_set(f,0,z-sin(x));
	return GSL_SUCCESS;
}

	double root(double z){
	const gsl_multiroot_fsolver_type* T =gsl_multiroot_fsolver_hybrids;
	gsl_multiroot_fsolver* S = gsl_multiroot_fsolver_alloc(T,1);

	gsl_multiroot_function F;
	F.f=root_equation; /* Struct F containing the function, its dimension and parameters */
	F.n=1;
	F.params=(void*)&z;

	gsl_vector* start = gsl_vector_alloc(1);
	gsl_vector_set(start,0,0);
	gsl_multiroot_fsolver_set(S,&F,start);



	int iter=0;
	int flag; /* precedure for checking of the steps are small */
	do{
	iter++;
	gsl_multiroot_fsolver_iterate(S);
	flag=gsl_multiroot_test_residual(S->f,1e-13);
	}

	while(flag==GSL_CONTINUE);

	double result = gsl_vector_get(S->x,0);
	gsl_vector_free(start);
	gsl_multiroot_fsolver_free(S);
	return result;
}


int main(){
        printf("i root(i) sqrt(i)\n");
        for(double i=-1;i<=1;i+=0.05)
        printf("%8.3f %8.3f %8.3f\n",i,root(i),asin(i));
}
