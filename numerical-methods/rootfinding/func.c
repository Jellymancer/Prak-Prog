#include "func.h"

void fRos(gsl_vector* p,gsl_vector* fx){
                double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1);
                gsl_vector_set(fx,0, 400*x*x*x-2*x*(200*y-1)-2);
                gsl_vector_set(fx,1, 100*2*(y-x*x));
}

void Jacobi_fRos(gsl_vector* p, gsl_matrix* J){
                double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1);
                gsl_matrix_set(J,0,0,2+100*2*(-2*x)*(-2*x)+100*2*(y-x*x)*(-1)*2);
                gsl_matrix_set(J,0,1,100*2*(-2*x));
                gsl_matrix_set(J,1,0,100*2*(-1)*2*x); gsl_matrix_set(J,1,1,100*2);}


void fSys(gsl_vector* p, gsl_vector* fx){
                double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1); double A=10000;
                gsl_vector_set(fx,0,A*x*y-1);
                gsl_vector_set(fx,1,exp(-x)+exp(-y)-1-1/A);
                }

void Jacobi_fSys(gsl_vector* p, gsl_matrix* J){
                double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1); double A=10000;
		gsl_matrix_set(J,0,0,A*y); gsl_matrix_set(J,0,1,A*x);
		gsl_matrix_set(J,1,0,-1*exp(-1*x)); gsl_matrix_set(J,1,1,-1*exp(-1*y));
}

void fHim(gsl_vector* p, gsl_vector* fx){
                double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1);
                gsl_vector_set(fx,0,4*x*x*x+2*x*(2*y-21)+2*(y*y-7));
                gsl_vector_set(fx,1,4*y*y*y+2*y*(2*x-13)+2*x*x-22);
                }


void Jacobi_fHim(gsl_vector* p, gsl_matrix* J){
                double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1);
                gsl_matrix_set(J,0,0,4*3*x*x+2*(2*y-21)); gsl_matrix_set(J,0,1,2*y*2+2*2*x);
                gsl_matrix_set(J,1,0,2*x*2+2*2*y); gsl_matrix_set(J,1,1,4*3*y*y+2*(2*x-13));
}

int GSLfRos(const gsl_vector* p, void *params,gsl_vector* f){
	double x = gsl_vector_get (p,0);
	double y = gsl_vector_get (p,1);
	double y0 = 1*(1 - x);
	double y1 = 100 * (y-x*x);

	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);

	return GSL_SUCCESS;
}

int GSLfSys(const gsl_vector* p, void *params,gsl_vector* f){
	double A=10000;
        double x = gsl_vector_get (p,0);
        double y = gsl_vector_get (p,1);
        double y0 = A*x*y-1;
        double y1 = exp(-x)+exp(-y)-(1+1/A);

        gsl_vector_set (f, 1, y0);
        gsl_vector_set (f, 0, y1);

        return GSL_SUCCESS;
}

int GSLfHim(const gsl_vector* p, void *params,gsl_vector* f){
        double x = gsl_vector_get (p,0);
        double y = gsl_vector_get (p,1);
        double y0 = x*x+y-11;
        double y1 = x+y*y-7;

        gsl_vector_set (f, 1, y0);
        gsl_vector_set (f, 0, y1);

        return GSL_SUCCESS;
}

