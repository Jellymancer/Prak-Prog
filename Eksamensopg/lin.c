#include<math.h>
#include<assert.h>
#include "interpol.h"

double linterp(int n, double* x, double* y, double z){
int i;
i=binsearch(n,z,x);
assert(n>1 && z>=x[0] && z<=x[n-1]); // assers if z is between the start and end of x.
        for(int j=0;j<n-1;j++){
                double xlow=x[j], xhigh=x[j+1];
                if(z<=xhigh & z>=xlow) i=j;
        }

	double ai=y[i]; // calculate the ai and bi for the line
	double bi=(y[i+1]-y[i])/(x[i+1]-x[i]);
	return ai+bi*(z-x[i]); //return the y value for the spline at point z.
}


double linterpint(int n, double* x, double* y, double z){
	int i;
	i=binsearch(n,z,x);

	double ai, bi, ak, bk;
	double intsum = 0.0;
	for(int k=0;k<=i-1;k++){ // find the integral of the spline between all points. The integral is the sum of the individual line
	//integrals between the points listed in x.
	ak=y[k];
	bk=(y[k+1]-y[k])/(x[k+1]-x[k]);
	intsum = intsum + ak*x[k+1]+1/2*pow(x[k+1],2)*bk - ak*x[k]-1/2*pow(x[k],2)*bk; //summing the integrals
	}

	ai=y[i]; //last integral evaluated seperately. 
	bi=(y[i+1]-y[i])/(x[i+1]-x[i]);
	intsum = intsum + ai*z+1/2*pow(z,2)*bi - ai*x[i]-1/2*pow(x[i],2)*bi; 
	return intsum;
}
