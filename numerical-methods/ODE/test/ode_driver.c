#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define float double
#define ode_stepper rkstep12
void ode_stepper(
        void f(int n,float x,float*y,float*dydx),
        int n, float x, float* y, float h, float* yh, float* dy);

int ode_driver(void f(int n,float x,float*y,float*dydx),
int n,float*xlist,float**ylist,
float b,float h,float acc,float eps,int max){
  int i,k=0; float x,*y,s,err,normy,tol,a=xlist[0],yh[n],dy[n];
  while(xlist[k]<b){
    x=xlist[k], y=ylist[k]; if(x+h>b) h=b-x;
    ode_stepper(f,n,x,y,h,yh,dy);
    s=0; for(i=0;i<n;i++) s+=dy[i]*dy[i]; err  =sqrt(s);
    s=0; for(i=0;i<n;i++) s+=yh[i]*yh[i]; normy=sqrt(s);
    tol=(normy*eps+acc)*sqrt(h/(b-a));
    if(err<tol){ /* accept step and continue */
      k++; if(k>max-1) return -k; /* uups */
      xlist[k]=x+h; for(i=0;i<n;i++)ylist[k][i]=yh[i];
      }
    if(err>0) h*=pow(tol/err,0.25)*0.95; else h*=2;
    } /* end while */
  return k+1; } /* return the number of entries in xlist/ylist */

