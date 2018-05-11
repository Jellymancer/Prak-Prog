#include "jac.h"

double* jac_onevalue(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int row,int ascending){

int changed, sweeps=0, n=A->size1;

for(int i=0;i<n;i++)
	gsl_vector_set(e,i,gsl_matrix_get(A,i,i));


gsl_matrix_set_identity(V);
double phi;

do{ changed=0; int p,q; sweeps++;
	p=row-1; for(q=p+1;q<n;q++){
		double app=gsl_vector_get(e,p);
		double aqq=gsl_vector_get(e,q);
		double apq=gsl_matrix_get(A,p,q);
		if(ascending==1)
			phi=0.5*atan2(2*apq,aqq-app);
		if(ascending==0)
			phi=-0.5*atan2(2*apq,app-aqq);

		double c = cos(phi), s = sin(phi);
		double app1=c*c*app-2*s*c*apq+s*s*aqq;
		double aqq1=s*s*app+2*s*c*apq+c*c*aqq;
		if(app1!=app || aqq1!=aqq){ changed=1;
			gsl_vector_set(e,p,app1);
			gsl_vector_set(e,q,aqq1);
			gsl_matrix_set(A,p,q,0.0);
			for(int i=0;i<p;i++){
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
			} } }while(changed!=0);

static double result[2];
result[0]=gsl_vector_get(e,row-1);
result[1]=sweeps;
return result;}


double jac_eigvlbyeigvl(gsl_matrix* A, gsl_vector* e, gsl_matrix* V){
int changed, sweeps=0, n=A->size1;
gsl_matrix* Aoriginal=gsl_matrix_calloc(n,n);
gsl_matrix_memcpy(Aoriginal,A);
double phi;
double totalsweeps=0;
gsl_matrix_set_identity(V);
gsl_matrix_memcpy(A,Aoriginal);
	for(int i=0;i<n;i++){
		gsl_vector_set(e,i,gsl_matrix_get(A,i,i));
	}
for(int p=0;p<n;p++){

do{ changed=0; int q; sweeps++;
			for(q=p+1;q<n;q++){
				double app=gsl_vector_get(e,p);
				double aqq=gsl_vector_get(e,q);
				double apq=gsl_matrix_get(A,p,q);
				phi=0.5*atan2(2*apq,aqq-app);
				double c = cos(phi), s = sin(phi);
				double app1=c*c*app-2*s*c*apq+s*s*aqq;
				double aqq1=s*s*app+2*s*c*apq+c*c*aqq;
			if(app1!=app || aqq1!=aqq){ changed=1;
			gsl_vector_set(e,p,app1);
			gsl_vector_set(e,q,aqq1);
			gsl_matrix_set(A,p,q,0.0);
			for(int i=0;i<p;i++){
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
			} } }while(changed!=0);



totalsweeps = totalsweeps+sweeps;
}
gsl_matrix_free(Aoriginal);
return totalsweeps;}
