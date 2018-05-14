#include "func.h"
#define MODEL(t) A*exp(-(t)/T)+B

void qrgsdecomp(gsl_matrix *E,gsl_matrix *W){
int s=E->size2; //Number of columns in matrix E
for(int i=0;i<s;i++){
        gsl_vector_view col=gsl_matrix_column(E,i);
        double Rii = gsl_blas_dnrm2(&col.vector);
        gsl_matrix_set(W,i,i,Rii);
        gsl_vector_scale(&col.vector,1/Rii);

        for(int j=i+1;j<s;j++){
                gsl_vector_view col2= gsl_matrix_column(E,j);
                double Rij = 0;
                gsl_blas_ddot(&col.vector,&col2.vector,&Rij);
                gsl_blas_daxpy(-Rij,&col.vector,&col2.vector);
                gsl_matrix_set(W,i,j,Rij);
                }}}

void qrgssolve(gsl_matrix* Q,  gsl_matrix* R,gsl_vector* b,gsl_vector* c){
gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,c);
for(int i=c->size-1; i>=0; i--){
        double s=gsl_vector_get(c,i);
        for(int k=i+1;k< c->size; k++)
                s-=gsl_matrix_get(R,i,k)*gsl_vector_get(c,k);
                gsl_vector_set(c,i,s/gsl_matrix_get(R,i,i));}}




double fRos(gsl_vector* p,gsl_vector* fx){
                double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1), fval;

		gsl_vector_set(fx,0, 400*x*x*x-2*x*(200*y-1)-2);
		gsl_vector_set(fx,1, 100*2*(y-x*x));

		fval= (1-x)*(1-x)+100* (y-x*x)*(y-x*x);
		return fval;
}

void Jacobi_fRos(gsl_vector* p, gsl_matrix* J){
                double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1);
                gsl_matrix_set(J,0,0,2+100*2*(-2*x)*(-2*x)+100*2*(y-x*x)*(-1)*2);
                gsl_matrix_set(J,0,1,100*2*(-2*x));
                gsl_matrix_set(J,1,0,100*2*(-1)*2*x); gsl_matrix_set(J,1,1,100*2);}



double fHim(gsl_vector* p, gsl_vector* fx){
                double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1), fval;
                gsl_vector_set(fx,0,4*x*x*x+2*x*(2*y-21)+2*(y*y-7));
                gsl_vector_set(fx,1,4*y*y*y+2*y*(2*x-13)+2*x*x-22);

		fval= pow(pow(x,2)+y-11,2) + pow(x+pow(y,2)-7,2);
		return fval;
                }


void Jacobi_fHim(gsl_vector* p, gsl_matrix* J){
                double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1);
                gsl_matrix_set(J,0,0,4*3*x*x+2*(2*y-21)); gsl_matrix_set(J,0,1,2*y*2+2*2*x);
                gsl_matrix_set(J,1,0,2*x*2+2*2*y); gsl_matrix_set(J,1,1,4*3*y*y+2*(2*x-13));

}

double fRad(gsl_vector* p, gsl_vector* fx){
                double A=gsl_vector_get(p,0), T=gsl_vector_get(p,1), B=gsl_vector_get(p,2), fval;
		double fxA=0,fxB=0,fxT=0,q,sum=0;

/*		double t[]= {0.47,1.41,2.36,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
		double y[]= {5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
		double e[]= {0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};
		int n = sizeof(t)/sizeof(t[0]);
*/

	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int n = sizeof(t)/sizeof(t[0]);

		for(int i=0;i<n;i++){ // setting derivatives and function to minimize F(A,T,B)=∑i(f(ti)-yi)²/σi²
			sum+=pow( (MODEL(t[i]) - y[i]) /e[i] ,2);
			q = 2*(MODEL(t[i])-y[i])/pow(e[i],2); // shorthand
			fxA += exp(-t[i]/T)*q;
			fxT += A*t[i]/pow(T,2)*exp(-t[i]/T)*q;
			fxB += q;
		}

		gsl_vector_set(fx,0,fxA);
		gsl_vector_set(fx,1,fxT);
		gsl_vector_set(fx,2,fxB);
                return sum;
                }

