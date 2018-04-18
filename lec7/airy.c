#include<stdio.h>
#include<gsl/gsl_sf_airy.h>
#include<math.h>


int main(){

	for(double x=-2*M_PI;x<0.2*M_PI;x+=0.013){
		double y = gsl_sf_airy_Ai(x,1);
		double y2 =  gsl_sf_airy_Bi(x,1);
		printf("%g %g %g\n",x,y,y2);
}

return 0;
}
