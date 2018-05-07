#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define float double // :)

int ode_driver(
        void f(int n, float x, float*y, float*dydx),
        int n, float* xlist, float** ylist,
        float b, float h, float acc, float eps, int max);

void f(int n, float x, float* y, float* dydx){
	dydx[0]=y[1];
	dydx[1]=-y[0];
	return;}

int main(){
	int n=2;
	int max=1000;
	float*xlist=(float*)calloc(max,sizeof(float));
	float**ylist=(float**)calloc(max,sizeof(float*)); 
	for(int i=0;i<max;i++) ylist[i]=(float*)calloc(n,sizeof(float)); 
	float pi=atan(1.)*4;
	float a=0, b=8*pi, h=0.1, acc=0.01, eps=0.01;
	xlist[0]=a; ylist[0][0]=0; ylist[0][1]=1;
	int k = ode_driver(f,n,xlist,ylist,b,h,acc,eps,max);
	if(k<0)printf("max steps reached in ode_driver\n");
	printf("# m=0, S=4\n");
	for(int i=0;i<k;i++)printf("%g %g\n",xlist[i],ylist[i][0]);
	printf("# m=1, S=0\n");
	for(int i=0;i<k;i++)printf("%g %g\n",xlist[i],sin(xlist[i]));
	return 0;
}
