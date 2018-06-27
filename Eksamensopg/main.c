#include "interpol.h"
#include<math.h>

int main(){
	int size = 10;
	/*define x,y points. Manual, as I want to showcase the wiggle reduction of akima splines
	and osplines.*/
	double x[10]={1,2,3,4,5,6,7,8,9,10};
	double y[10]={1,2,3,2,3,5,7,3,1,10};
	double df[10]={1,1,-1,1,2,2,-4,-2,9};
	//another x,y,p set is used, here p_i (df) is estimated using a quadratic spline. 
	double x3[10]={1,2,3,4,5,6,7,8,9,10};
	double y3[10]={1,30,1,1,1,1,1,-30,1,1};
	double df3[10];
	//allocatng splines.
	aspline* Q=aspline_alloc(size,x,y);
	cspline* Qc=cspline_alloc(size,x,y);
	aspline* Q3=aspline_alloc(size,x3,y3);
	cspline* Qc3=cspline_alloc(size,x3,y3);
	qspline* Qq3=qspline_alloc(size,x3,y3);

	FILE* manual;
	manual=fopen("manual.dat","w");
	for(int i=0;i<size;i++) fprintf(manual,"%g %g\n",x[i],y[i]);
	fprintf(manual,"\n\n");
	for(int i=0;i<size;i++){fprintf(manual,"%g %g\n",x3[i],y3[i]);
		df3[i]=qspline_deriv(Qq3,x3[i]);}

	ospline* Qo=ospline_alloc(size,x,y,df);
	ospline* Qo3=ospline_alloc(size,x3,y3,df3);

	double dz=0.01;
	fprintf(manual,"\n\n");
	for(double z=x[0];z<=x[size-1];z+=dz){
		double l=aspline_eval(Q,z);
		double c=cspline_eval(Qc,z);
		double o=ospline_eval(Qo,z);
		fprintf(manual,"%g %g %g %g\n",z,l,c,o);}
	fprintf(manual,"\n\n");
	for(double z=x3[0];z<=x3[size-1];z+=dz){
		double l=aspline_eval(Q3,z);
		double c=cspline_eval(Qc3,z);
		double o=ospline_eval(Qo3,z);
		fprintf(manual,"%g %g %g %g\n",z,l,c,o);}
	fclose(manual);

	//Now I make a sine curve in order to showcase the derivatives and integrals of akima.
	int size2=8;
	double x2[size2];
	double y2[size2];
	double df2[size2];
	for(int i=0;i<size2;i++){
		x2[i]=i;
		y2[i]=sin(i);
		df2[i]=cos(i);}

	qspline* Qq2=qspline_alloc(size2,x2,y2);

	aspline* Q2=aspline_alloc(size2,x2,y2);
	cspline* Qc2=cspline_alloc(size2,x2,y2);
	ospline* Qo2=ospline_alloc(size2,x2,y2,df2);

	FILE* sinsplines;
	sinsplines=fopen("sinspline.dat","w");
	for(double z=x2[0];z<=x2[size2-1];z+=dz){
		double c=cspline_eval(Qc2,z);
		double l=aspline_eval(Q2,z);
		double o=ospline_eval(Qo2,z);
		double dc=cspline_deriv(Qc2,z);
		double dl=aspline_deriv(Q2,z);
		double dog=ospline_deriv(Qo2,z);
		double d2c=cspline_deriv2(Qc2,z);
		double d2l=aspline_deriv2(Q2,z);
		double d2o=ospline_deriv2(Qo2,z);
		fprintf(sinsplines,"%g %g %g %g %g %g %g %g %g %g\n",z,c,l,o,dc,dl,dog,d2c,d2l,d2o);}
	fprintf(sinsplines,"\n\n");
	for(int i=0;i<size2;i++) fprintf(sinsplines,"%g %g\n",x2[i],sin(i));
	fprintf(sinsplines,"\n\n");
	for(int i=0;i<size2;i++) fprintf(sinsplines,"%g %g\n",x2[i],cos(i));
	fprintf(sinsplines,"\n\n");
	for(int i=0;i<size2;i++) fprintf(sinsplines,"%g %g\n",x2[i],-sin(i));
	fclose(sinsplines);


// integration of sine.

	double zmax=4;
	double cinteg = cspline_integ(Qc2,zmax);
	double ainteg = aspline_integ(Q2,zmax);
	double ointeg = ospline_integ(Qo2,zmax);
	printf("Integ cspline = %g\n Integ aspline =%g\n Integ cs-spline =%g\n",cinteg,ainteg,ointeg);
	return 0;
}

