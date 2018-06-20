#include "jac.h"
int jac(gsl_matrix* A, gsl_vector* e, gsl_matrix* V){

int changed, sweeps=0, n=A->size1;

for(int i=0;i<n;i++)
	gsl_vector_set(e,i,gsl_matrix_get(A,i,i));


gsl_matrix_set_identity(V);
do{ changed=0; int p,q;
	for(p=0;p<n;p++)for(q=p+1;q<n;q++){
		double app=gsl_vector_get(e,p);//finding the values to be rotated
		double aqq=gsl_vector_get(e,q);
		double apq=gsl_matrix_get(A,p,q);
		double phi=0.5*atan2(2*apq,aqq-app); //determining the rotation angle
		double c = cos(phi), s = sin(phi);
		double app1=c*c*app-2*s*c*apq+s*s*aqq;// used for comparing if values cahnge.
		double aqq1=s*s*app+2*s*c*apq+c*c*aqq;
		if(app1!=app || aqq1!=aqq){ changed=1; sweeps++; // checking if the eigenvalues change after rotation. if not, make next rotation.
			gsl_vector_set(e,p,app1); // setting new values for pp and qq indexes.
			gsl_vector_set(e,q,aqq1);
			gsl_matrix_set(A,p,q,0.0);
			for(int i=0;i<p;i++){ //performing the rotation for all involved numbers. method is described in the lecture book.
				double aip=gsl_matrix_get(A,i,p);
				double aiq=gsl_matrix_get(A,i,q);
				gsl_matrix_set(A,i,p,c*aip-s*aiq);
				gsl_matrix_set(A,i,q,c*aiq+s*aip); }
			for(int i=p+1;i<q;i++){
				double api=gsl_matrix_get(A,p,i);
				double aiq=gsl_matrix_get(A,i,q);
				gsl_matrix_set(A,p,i,c*api-s*aiq);
				gsl_matrix_set(A,i,q,c*aiq+s*api); }
			for(int i=q+1;i<n;i++){
				double api=gsl_matrix_get(A,p,i);
				double aqi=gsl_matrix_get(A,q,i);
				gsl_matrix_set(A,p,i,c*api-s*aqi);
				gsl_matrix_set(A,q,i,c*aqi+s*api); }
			for(int i=0;i<n;i++){
				double vip=gsl_matrix_get(V,i,p);
				double viq=gsl_matrix_get(V,i,q);
				gsl_matrix_set(V,i,p,c*vip-s*viq);
				gsl_matrix_set(V,i,q,c*viq+s*vip); }
			} } }while(changed!=0); // stop when rotation does not change eigenvalue.
return sweeps; }
