#include<math.h>
#include<assert.h>
#include "interpol.h"

double linterp(int n, double* x, double* y, double z){
	assert(n>1 && z>=x[0] && z<=x[n-1]);
	int i=0, j=n-1; /* binary search: */
	while(j-i>1){
		int m=(i+j)/2;
		if(z>=x[m]) i=m;
		else j=m;
	}
	double ai=y[i];
	double bi=(y[i+1]-y[i])/(x[i+1]-x[i]);
	return ai+bi*(z-x[i]);
}


double linterpint(int n, double* x, double* y, double z){
        assert(n>1 && z>=x[0] && z<=x[n-1]);
        int i=0, j=n-1; /* binary search: */

        while(j-i>1){
                int m=(i+j)/2;
                if(z>=x[m]) i=m;
                else j=m;
        }

	double ai, bi, ak, bk;
	double intsum = 0.0;
	for(int k=0;k<=i-1;k++){
	ak=y[k];
	bk=(y[k+1]-y[k])/(x[k+1]-x[k]);
	intsum = intsum + ak*x[k+1]+1/2*pow(x[k+1],2)*bk - ak*x[k]-1/2*pow(x[k],2)*bk;
	}

	ai=y[i];
	bi=(y[i+1]-y[i])/(x[i+1]-x[i]);
	intsum = intsum + ai*z+1/2*pow(z,2)*bi - ai*x[i]-1/2*pow(x[i],2)*bi; 
	return intsum; 
}
