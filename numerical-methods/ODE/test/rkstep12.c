void rkstep12(void f(int n, double x, double*yx, double*dydx),
int n, double x, double* yx, double h, double* yh, double* dy){
  double k0[n],yt[n],k12[n]; /* VLA: gcc -std=c99 */
  f(n,x    ,yx,k0);  for(int i=0;i<n;i++) yt[i]=yx[i]+ k0[i]*h/2;
  f(n,x+h/2,yt,k12); for(int i=0;i<n;i++) yh[i]=yx[i]+k12[i]*h;
  for(int i=0;i<n;i++) dy[i]=(k0[i]-k12[i])*h/2;
}
