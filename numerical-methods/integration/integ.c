#include "integ.h"


double integaux(
double f(double,double),  /*integration function*/
double a, /* start of integral*/
double b,  /*end of integral*/
double c,  /*used for infinite integration*/
double abstol,  /*absolute tolerance*/
double reltol,  /*relativ tolerance */
double f2,
double f3, /*subdivision points*/
int *nrecs) /*no of recursions*/

{
	double f1 = f(a+(b-a)/6,c);
	double f4 = f(a+5*(b-a)/6,c);
	double Q_high, Q_low; //Quadratures for low and high order
	Q_high=(2*f1+f2+f3+2*f4)/6*(b-a);
	Q_low=(f1+f4+f2+f3)/4*(b-a);

	double tau=abstol+reltol*fabs(Q_high);
	double err=fabs(Q_high-Q_low); // tolerance and error. Accept step if  err<tau...

	if(err < tau) return Q_high;

	else{/* Splits integral into two half integrals and recursively calls the intgrator again*/
		*nrecs = *nrecs+1;
		double subQ1=integaux(f,(a+b)/2,b,c,abstol,reltol,f1,f2,nrecs);
		double subQ2=integaux(f,a,(b+a)/2,c,abstol,reltol,f3,f4,nrecs);
	double newQ=subQ1+subQ2;
	return newQ;
	}
}


double integmain(double f(double,double), double a, double b, double abstol, double reltol, int *nrecs){
	double Q, f3,f2,a2,b2,c;
	double inf1(double x,double c){return f(x/(1-pow(x,2)),c) * ((1+pow(x,2))/pow(1-pow(x,2),2));}
	double inf2(double x, double c){return f(c+(x/(1-x)),c)*(1/(pow(1-x,2)));}
	double inf3(double x,double c){return f(c+(x/(1-x)),c)*(1/(pow(1+x,2)));}

	if(isinf(a)==-1 && isinf(b)==1){
		c=0;
		a2=-1;b2=1;
		f3 = inf1(a2+4*(b2-a2)/6,c);
		f2 = inf1(a2+2*(b2-a2)/6,c); /* weigths (for both Q_high and Q_low) and points are chosen
		so that the points  are reused when calcutaling the riemann integrals after consecutive subdivisions.
		The used weigths and points are taken from the ingration.pdf book for the course*/
		Q =  integaux(inf1, a2, b2,c, abstol, reltol, f2, f3, nrecs);
	}
	else if(isinf(a)==0 && isinf(b)==1){
		c=a;
		a2=0;b2=1;
		f3 = inf2(a2+4*(b2-a2)/6,c);
		f2 = inf2(a2+2*(b2-a2)/6,c);
		Q =  integaux(inf2, a2, b2,c, abstol, reltol, f2, f3, nrecs);
	}
	else if(isinf(a)==-1 && isinf(b)==0){
		c=b;
		a2=-1;b2=0;
		f3 = inf3(a2+4*(b2-a2)/6,c);
		f2 = inf3(a2+2*(b2-a2)/6,c);
		Q =  integaux(inf3, a2, b2,c, abstol, reltol, f2, f3, nrecs);
	}
	else{
		c=0;
		f3 = f(a+4*(b-a)/6,c);
		f2 = f(a+2*(b-a)/6,c);
		Q =  integaux(f, a, b,c, abstol, reltol, f2, f3, nrecs);
	}

	return Q;
}
