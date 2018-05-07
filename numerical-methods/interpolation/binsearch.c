#include "interpol.h"

int binsearch(int n, double z, double* x){
//this function evaluates between which points (listed in x) z lies.
assert(n>1 && z>=x[0] && z<=x[n-1]); // assers if z is between the start and end of x.
	for(int i=0;i<n;i++){
		double xlow=x[i], xhigh=x[i+1];
		if(z<=xhigh & z>=xlow) return i;
	}
}

int binsearchqspl(int n, double z, qspline* s){
//this function evaluates between which points (listed in x) z lies.
assert(n>1 && z>=s->x[0] && z<=s->x[n-1]); // assers if z is between the start and end of x.
	for(int i=0;i<n;i++){
		double xlow=s->x[i], xhigh=s->x[i+1];
		if(z<=xhigh & z>=xlow) return i;
	}
}


