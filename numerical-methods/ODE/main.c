#include "print.h"
#include "functions.h"


void f_trial(double x, gsl_vector* y, gsl_vector* dydx){
	gsl_vector_set(dydx,0,3*x+4);
}

void f_airy(double x, gsl_vector* y, gsl_vector* dydx){
	double y0, y1;
	y0=gsl_vector_get(y,0);
	y1=gsl_vector_get(y,1);
	gsl_vector_set(dydx,0,y1);
	gsl_vector_set(dydx,1,y0*x);}

int main(){
// trial function dydx = 3x+4
printf("Assignment 1: ODE using Runge-Kutta midpoint method\nTesting the algorythm on a simple differential equation: dy/dx = 3x+4 which has an analytic solution y= 3/2 x² + 4x + c.\n");
int n=1; // size. Has to be equal to number of equations to solve.
gsl_vector* yx=gsl_vector_alloc(n);

double x_start=-50; // startpoint x
double step = 0.1; //stepsize
double y0start=3*(50*50)*0.5-4*50; // start values for the dy/dx.
int max=10000; //max no of iterations.
gsl_vector_set(yx,0,y0start);
printf("Setting startpoint x=-50 and dydx=3/2*50²-200 so that c=0.\n");
double abstol=0.001, reltol=0.001, x_end=50; //  absolute tolerance, relative tolerance and end points


printf("I also set the end point to be x_end = 50. The analytical solution at this point gives y(50)=3950.\n");
int iter=driver(x_start,x_end,step,yx,abstol,reltol,max,rkstepX,f_trial);
printf("The obtained numerical solution is: y(50)=%f.\n",gsl_vector_get(yx,0));
printf("The number of iterations made is: %d\n",iter);

printf("\n\nAssignment 2: Storing the path.\n Using the same function, the algorythm is performed again, this time the points calcutated along the way to x=50 are stored.\n");
printf("The numerical and analytical solutions are shown in figure 1");
printf("\nThe numerical solution alone is plotted in figure 2 using a lower initial step size and tolerances. It can be seem that steps become denser close to 0 where the derivative changes rapidly.");
gsl_matrix* xypath=gsl_matrix_alloc(max,n+1);
gsl_matrix* xypathL=gsl_matrix_alloc(max,n+1);

x_start=-50;  gsl_vector_set(yx,0,y0start);
iter = driverwpath(x_start,x_end,step*100,yx,abstol,reltol,max,rkstepX,f_trial,xypath,1);
FILE* as1; as1=fopen("as1.txt","w");
gsl_matrix_view xycutpath=gsl_matrix_submatrix(xypath,0,0,iter,n+1);
printm(&xycutpath.matrix,as1);
fprintf(as1,"\n\n");

gsl_vector_set(yx,0,y0start);
x_start=-50;  gsl_vector_set(yx,0,y0start);
iter = driverwpath(x_start,x_end,step,yx,abstol*100,reltol*100,max,rkstepX,f_trial,xypathL,1);
gsl_matrix_view xycutpathL=gsl_matrix_submatrix(xypathL,0,0,iter,n+1);
printm(&xycutpathL.matrix,as1);
fclose(as1);
FILE* as1ana; as1ana=fopen("as1ana.txt","w"); // analytical solution to the teste diff equation is printed
for(double w=x_start ;w<=x_end; w+=2){
		double  anaf = 3*pow(w,2)*0.5+4*w;
		fprintf(as1ana,"%g %g\n",w,anaf);
}
fclose(as1ana);



// Airy functions
printf("\n\nNow I solve the differential equation d²y/dx² -xy = 0. The solution to this equation are the airy functions.\n");
printf("The equation is split into dy1/dx=y2 and dy2/dx =xy1\nTwo functions solve this equations the Ai and Bi airy functions.\n");
printf("I use the Bi starting values (from wikipedia) so the program should find Bi. The numeric solution is plotted together with the Bi from GSL in figure 3\n");

n=2;
gsl_vector* yxairy=gsl_vector_alloc(n);
gsl_matrix* xypath2pos=gsl_matrix_alloc(max,n+1);
gsl_matrix* xypath2neg=gsl_matrix_alloc(max,n+1);

FILE* GSLa2; GSLa2=fopen("GSLairy.txt","w"); // prints analytical Airy function.
for(double w=-10 ;w<=3; w+=0.5){
		double airy1 = gsl_sf_airy_Ai(w,1);
		double airy2 =  gsl_sf_airy_Bi(w,1);
		fprintf(GSLa2,"%g %g %g\n",w,airy1,airy2);
}
fclose(GSLa2);
//Airy functions are split into the negative and positive parts, since I only know the limiting
//values at x=0.
x_start=0;
double y1_start=0.6149266, y2_start=0.4482883;
gsl_vector_set(yxairy,0,y1_start);
gsl_vector_set(yxairy,1,y2_start);
max=10000;
abstol=0.001; reltol=0.001; //  absolute tolerance, relative tolerance and end points

//Positive first.
step=0.1; x_end=3;
int iterpos = driverwpath(x_start,x_end,step,yxairy,abstol,reltol,max,rkstepX,f_airy,xypath2pos,1);
//and negative (taken in negative direction from 0 to -10);
gsl_vector_set(yxairy,0,y1_start);
gsl_vector_set(yxairy,1,y2_start);
step=-0.1; x_end=-10;
int iterneg = driverwpath(x_start,x_end,step,yxairy,abstol,reltol,max,rkstepX,f_airy,xypath2neg,0);

//The solutions are combined and printed
gsl_matrix_view xycutpath2pos=gsl_matrix_submatrix(xypath2pos,0,0,iterpos,n+1);
gsl_matrix_view xycutpath2neg=gsl_matrix_submatrix(xypath2neg,0,0,iterneg,n+1); //remove all 0-values (matrix is larger than needed)


FILE* as2; as2=fopen("airy.txt","w");
printm(&xycutpath2neg.matrix,as2);
printm(&xycutpath2pos.matrix,as2);
fclose(as2);


return 0;
}
