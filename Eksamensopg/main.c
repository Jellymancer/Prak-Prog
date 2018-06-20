#include "interpol.h"
#include<math.h>

int main(){
	int size = 15;
	/*define x,y points. Manual, as I want to showcase the wiggle reduction of akima splines
	and osplines.*/
	double x[15]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
	double y[15]={1,1,1,1,30,1,1,1,1,1,8,8,8,8,8};
	double df[15]; //allocate x, y and df (derivatives only used for osplines) points is lists.
	aspline* Q=aspline_alloc(size,x,y);
	cspline* Qc=cspline_alloc(size,x,y);
	qspline* Qq=qspline_alloc(size,x,y);

	for(int i=0;i<size;i++){ printf("%g %g\n",x[i],y[i]);
		df[i]=cspline_deriv(Qc,x[i]);
		printf("%g\n",df[i]);
	}
	ospline* Qo=ospline_alloc(size,x,y,df);

	FILE* asplines;
	asplines=fopen("aspline.dat","w");
	double dz=0.01;
	for(double z=x[0];z<=x[size-1];z+=dz){
		double l=aspline_eval(Q,z);
		double c=cspline_eval(Qc,z);
		double o=ospline_eval(Qo,z);
		fprintf(asplines,"%g %g %g %g\n",z,l,c,o);}
	fclose(asplines);

	//Now I make a sine curve in order to showcase the derivatives and integrals of akima.
	int size2=15;
	double x2[size2];
	double y2[size2];
	double df2[size2];

	for(int i=1;i<size2;i++){
		x2[i]=i;
		y2[i]=sin(0.5*i);}
	qspline* Qq2=qspline_alloc(size2,x2,y2);
	for(int i=0;i<size2-1;i++){df2[i]=qspline_deriv(Qq2,x2[i]);printf("%g\n",df2[i]);}

	aspline* Q2=aspline_alloc(size2,x2,y2);
	cspline* Qc2=cspline_alloc(size2,x2,y2);
	ospline* Qo2=ospline_alloc(size2,x2,y2,df2);

	FILE* sinsplines;
	sinsplines=fopen("sinspline.dat","w");
	for(double z=x2[0];z<=x2[size2-1];z+=dz){
		double dc=cspline_deriv(Qc2,z);
		double dl=aspline_deriv(Q2,z);
		double dog=ospline_eval(Qo2,z);
		double d2c=cspline_deriv2(Qc2,z);
		double d2l=aspline_deriv2(Q2,z);
		fprintf(sinsplines,"%g %g %g %g %g %g %g\n",z,dc,dl,dog,d2c,d2l);}
	fprintf(sinsplines,"\n\n");
	for(int i=0;i<size2;i++) fprintf(sinsplines,"%g %g\n",x2[i],0.5*cos(0.5*i));
	fprintf(sinsplines,"\n\n");
	for(int i=0;i<size2;i++) fprintf(sinsplines,"%g %g\n",x2[i],-0.25*sin(0.5*i));
	fclose(sinsplines);
	double zmax=4;
	double cinteg = cspline_integ(Qc2,zmax);
	double ainteg = aspline_integ(Q2,zmax);
//	printf("Integ cspline = %g\n Integ aspline =%g\n",cinteg,ainteg); 
	return 0;
}

