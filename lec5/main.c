#include "nvector.h"
#include "stdio.h"
#include "stdlib.h"
#define RND (double)rand()/RAND_MAX


int main(){
	int n=4;

	printf("\nTesting _aloc, _set, _get and _print.\nAllocating a vector of size 4:\n");
	nvector *v = nvector_alloc(n);
	if(v==NULL) fprintf(stderr,"Failed to allocate vector");
	else{printf("Allocation sucessfull.\ All elements are now set to 1 and retrived.\n");
	printf("Checking if the values of the vector equal one after runnig _set by using _get.\n");
	for(int i=0; i<n;i++){ nvector_set(v,i,1);
	printf("Value at index i=%d is equal to %g\n",i,nvector_get(v,i));}
	printf("\nNow testing the print function which should return the same vector of ones.\n");
	nvector_print("v=\n",v);}

	printf("\n\nTesting _set_zero by allocating new vector and checking if all entries are zero.\n");
	nvector_set_zero(v);
	nvector_print("v=\n",v);


	printf("\n\nNow I define two vectors, a and b:\n");
	nvector *a = nvector_alloc(n);
	nvector *b = nvector_alloc(n);
	nvector *c = nvector_alloc(n);
	for(int i=0; i<n;i++){ nvector_set(a,i,round(RND+RND+RND+RND+RND)); nvector_set(b,i,round(RND+RND+RND)); }
	nvector_print("a=\n",a);
	printf("\n");
	nvector_print("b=\n",b);
	printf("\nUsing _add and _sub to subtract and add the vectors yields:\n");

	nvector_add(a,b,c);
	nvector_print("a+b=\n",c);
	printf("\n");
	nvector_sub(a,b,c);
	nvector_print("a-b=\n",c);

	printf("\n\n Defining x=5 and checking b*x:\n");
	double x=5;
	nvector_scale(b,c,x);
	nvector_print("x*b=\n",c);

return 0;
}
