void rkstep23( void f(int n,double x,double* y,double* dydx),
int n, double x, double* yx, double h, double* yh, double* dy){
  int i; double k1[n],k2[n],k3[n],k4[n],yt[n]; /* VLA: -std=c99 */
  f(n,   x    ,yx,k1); for(i=0;i<n;i++) yt[i]=yx[i]+1./2*k1[i]*h;
  f(n,x+1./2*h,yt,k2); for(i=0;i<n;i++) yt[i]=yx[i]+3./4*k2[i]*h;
  f(n,x+3./4*h,yt,k3); for(i=0;i<n;i++)
    yh[i]=yx[i]+(2./9 *k1[i]+1./3*k2[i]+4./9*k3[i])*h;
  f(n,  x+h   ,yh,k4);  for(i=0;i<n;i++){
    yt[i]=yx[i]+(7./24*k1[i]+1./4*k2[i]+1./3*k3[i]+1./8*k4[i])*h;
    dy[i]=yh[i]-yt[i];
    }
}

