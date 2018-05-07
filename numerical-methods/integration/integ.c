#include<stdio.h>
#include<math.h>
#include<assert.h>

double integ(double f(double),  /*integration function*/
double a, /* start of integral*/
double b,  /*end of integral*/
double abstol,  /*absolute tolerance*/
double reltol,  /*relativ tolerance */
double *err /*error estimate*/)
{
	double Q_high, Q_low; //Quadratures for low and high order
	Q_high = (2*f1+f2+f3+2*f4)/6 * (b-a);
	Q_low= f4=f(a+5*(b-a)/6);

	tau=abstol+reltol*fabs(Q_high);
	err=fabs(Q-q); // tolerance and error. Accept step if  err<tau...
	if(error < tolerance) return Q; 
	else{/* Splits integral into two half integrals and recursively calls the intgrator again */
		double 	subQ1=integ(f,(a+b)/2,b,abstol,reltol,*err1)
		double 	subQ2=integ(f,a,(b+a)/2,abstol,reltol,*err2)
		*err=sqrt(power(err1,2)+power(err2,2)); //error found for the new integrals.
	double newQ=subQ1+subQ2;
	return newQ;
	}
}
