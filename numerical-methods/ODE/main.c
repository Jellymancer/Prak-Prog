#include "print.h"


void f_trial(double x, gsl_vector* y, gsl_vector* dydx){
	gsl_vector_set(dydx,0,3*x+4);
}

void f_airy(double x, gsl_vector* y, gsl_vector* dydx){
	double dydx0, dydx1;
	dydx0=gsl_vector_get(y,0);
	dydx1=gsl_vector_get(y,1);

	gsl_vector_set(dydx,0,dydx1);
	gsl_vector_set(dydx,1,x*dydx0);}

#include "functions.h"



int main(){
// trial function dydx = 3x+4
printf("Assignment 1: ODE using Runge-Kutta midpoint method\nTesting the algorythm on a simple differential equation: dy/dx = 3x+4 which has an analytic solution y= 3/2 x² + 4x + c.\n");
int n=1; // size. Has to be equal to number of equations to solve.
gsl_vector* yx=gsl_vector_alloc(n);

double x_start=0; // startpoint x
double step = 0.01; //stepsize
double y0start=0; // start values for the dy/dx.
int max=10000; //max no of iterations.
gsl_vector_set(yx,0,y0start);
printf("Setting startpoint x=0 and dydx=0 so that c=0.\n");
double abstol=0.001, reltol=0.001, x_end=50; //  absolute tolerance, relative tolerance and end points


printf("I also set the end point to be x_end = 50. The analytical solution at this point gives y(50)=3950.\n");
int iter=driver(x_start,x_end,step,yx,abstol,reltol,max,rkstepX,f_trial);
printf("The obtained numerical solution is: y(50)=%f.\n",gsl_vector_get(yx,0));
printf("\n%d\n",iter);

printf("\n\nAssignment 2: Storing the path.\n Using the same function, the algorythm is performed again, this time the points calcutated along the way to x=50 are stored.\n");
printf("The numerical and analytical solutions are shwon in figure 1");
gsl_matrix* xypath=gsl_matrix_alloc(max,n+1);
x_start=0;  gsl_vector_set(yx,0,y0start);
iter = driverwpath(x_start,x_end,step,yx,abstol,reltol,max,rkstepX,f_trial,xypath);
FILE* as1; as1=fopen("as1.txt","w");
gsl_matrix_view xycutpath=gsl_matrix_submatrix(xypath,0,0,iter,n+1);
printm(&xycutpath.matrix,as1);
fclose(as1);



// Airy functions
printf("\n\nNow I solve the differential equation d²y/dx² -xy = 0. The solution to this equation are the airy function.\n");
printf("The equation is split into dy1/dx=y2 and dy2/dx =xy1");
n=2;
gsl_vector* yxairy=gsl_vector_alloc(n);
gsl_matrix* xypath2=gsl_matrix_alloc(max,n+1);

x_start=-6.28319; x_end=1;
FILE* GSLa2; GSLa2=fopen("GSLairy.txt","w");
for(double w=x_start ;w<x_end; w+=0.2){
		double airy1 = gsl_sf_airy_Ai(w,1);
		double airy2 =  gsl_sf_airy_Bi(w,1);
		fprintf(GSLa2,"%g %g %g\n",w,airy1,airy2);
}
fclose(GSLa2);

double y1_start=gsl_sf_airy_Bi(x_start,1), y2_start=gsl_sf_airy_Ai(x_start,1);
gsl_vector_set(yxairy,0,y1_start);
gsl_vector_set(yxairy,1,y2_start);
max=10000; step = 0.0001;
abstol=0.001; reltol=0.001; //  absolute tolerance, relative tolerance and end points

iter = driverwpath(x_start,x_end,step,yxairy,abstol,reltol,max,rkstepX,f_airy,xypath2);

FILE* as2; as2=fopen("airy.txt","w");
gsl_matrix_view xycutpath2=gsl_matrix_submatrix(xypath2,0,0,iter,n+1);
printm(&xycutpath2.matrix,as2);
fclose(as2);


return 0;
}
