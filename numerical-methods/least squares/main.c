#include "matrixprep.h"
#include "fit.h"

void fittry(gsl_vector* c, double* x, double* y, double* dy, double m, gsl_vector* dc){
                double dx=0.1;
                double fx=0, fxplus=0, fxminus=0;
                for(double xx=x[0];xx<x[8];xx+=dx){
                for(int k=0;k<m;k++){
                        fx+=funs(k,xx)*gsl_vector_get(c,k);
                	fxplus+=funs(k,xx)*gsl_vector_get(c,k) + gsl_vector_get(dc,k);
                	fxminus+=funs(k,xx)*gsl_vector_get(c,k) - gsl_vector_get(dc,k);
		}
                       printf("%f %f %f %f\n",xx,fx,fxminus,fxplus);
                	fx=0.0; fxminus=0.0; fxplus=0.0;}
}
int main(){
double x[]  ={0.1  ,  1.33,    2.55  ,  3.78   ,    5  ,  6.22  ,  7.45  ,  8.68  ,   9.9};
double y[]  = {-15.3,    0.32  ,  2.45   , 2.75    ,2.27   , 1.35 ,  0.157  , -1.23 ,  -2.75};
double dy[] = {1.04,   0.594  , 0.983   ,0.998    ,1.11   ,0.398  , 0.535  , 0.968 ,  0.478};
int n=sizeof(x)/sizeof(x[0]);
FILE* ogdata;
ogdata=fopen("ogdata.txt","w");
for(int i=0;i<n;i++) fprintf(ogdata,"%f %f %f\n",x[i],y[i],dy[i]);
fclose(ogdata);
int m=3;
gsl_matrix* M=gsl_matrix_calloc(n,m);
gsl_vector* b=gsl_vector_calloc(n);
gsl_vector* c=gsl_vector_calloc(m);
gsl_matrix* Sigma=gsl_matrix_calloc(m,m);
matrix_generation(M,b,x,y,dy);
lsfit(M,b,c,Sigma);

gsl_vector* dc = gsl_vector_alloc(m);
	for(int k=0;k<m;k++){
		double skk=gsl_matrix_get(Sigma,k,k);
		gsl_vector_set(dc,k,sqrt(skk));
		}
printm(Sigma,stdout);
printf("\n\n\n\n\n");
fittry(c,x,y,dy,m,dc);
FILE* opg2; opg2 = fopen("opg2.txt","w");
fprintf(opg2,"The coefficients for the fit: c_0*log(x)+c_1+c_2*x found via least squares method are\n");
printv(c,opg2);
fprintf(opg2,"\nThe covariance matrix for the coefficients is calculated to be:\n");
printm(Sigma,opg2);
fprintf(opg2,"Best fit is plotted together with the maximum/minimum fits (using maximum/minimum deviations of the fit parameters)\n");
return 0;}



