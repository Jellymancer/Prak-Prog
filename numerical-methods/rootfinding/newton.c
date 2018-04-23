#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>

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
                }
        }

}


void qrgssolve(gsl_matrix* Q,  gsl_matrix* R,gsl_vector* b,gsl_vector* c){
gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,c);
for(int i=c->size-1; i>=0; i--){
        double s=gsl_vector_get(c,i);
        for(int k=i+1;k< c->size; k++)
                s-=gsl_matrix_get(R,i,k)*gsl_vector_get(c,k);
                gsl_vector_set(c,i,s/gsl_matrix_get(R,i,i));}

}




void newton(void f(gsl_vector* x,gsl_vector* fx, gsl_matrix *J), gsl_vector* x, double dx, double eps){
	int n=x->size;
	gsl_matrix* J = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* z  = gsl_vector_alloc(n);
	gsl_vector* fz = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);
	while(1){
		f(x,fx,J);
		for (int j=0;j<n;j++){
			gsl_vector_set(x,j,gsl_vector_get(x,j)+dx);
			f(x,df,J);
			gsl_vector_sub(df,fx); /* df=f(x+dx)-f(x) */
			for(int i=0;i<n;i++) gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/dx);
			gsl_vector_set(x,j,gsl_vector_get(x,j)-dx);
			}
		qrgsdecomp(J,R);
		qrgssolve(J,R,fx,Dx);
		gsl_vector_scale(Dx,-1);
		double s=1;
		while(1){
			gsl_vector_memcpy(z,x);
			gsl_vector_add(z,Dx);
			f(z,fz,J);
			if( gsl_blas_dnrm2(fz)<(1-s/2)*gsl_blas_dnrm2(fx) || s<0.02 ) break;
			s*=0.5;
			gsl_vector_scale(Dx,0.5);
			}
		gsl_vector_memcpy(x,z);
		gsl_vector_memcpy(fx,fz);
		if( gsl_blas_dnrm2(Dx)<dx || gsl_blas_dnrm2(fx)<eps ) break;
		}
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(fx);
	gsl_vector_free(z);
	gsl_vector_free(fz);
	gsl_vector_free(df);
	gsl_vector_free(Dx);
}

void newtonana(void f(gsl_vector* x,gsl_vector* fx,gsl_matrix* J), gsl_vector* x, double dx, double eps,gsl_matrix* J){
	int n=x->size;
	gsl_matrix* Q = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* z  = gsl_vector_alloc(n);
	gsl_vector* fz = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);

	while(1){
		f(x,fx,J);
		for (int j=0;j<n;j++){
			gsl_vector_set(x,j,gsl_vector_get(x,j)+dx);
			f(x,df,J);
			gsl_vector_sub(df,fx); /* df=f(x+dx)-f(x) */
			gsl_vector_set(x,j,gsl_vector_get(x,j)-dx);
			}
		gsl_matrix_memcpy(Q,J);
		qrgsdecomp(Q,R);
		qrgssolve(Q,R,fx,Dx);
		gsl_vector_scale(Dx,-1);
		double s=1;
		while(1){
			gsl_vector_memcpy(z,x);
			gsl_vector_add(z,Dx);
			f(z,fz,J);
			if( gsl_blas_dnrm2(fz)<(1-s/2)*gsl_blas_dnrm2(fx) || s<0.02 ) break;
			s*=0.5;
			gsl_vector_scale(Dx,0.5);
			}
		gsl_vector_memcpy(x,z);
		gsl_vector_memcpy(fx,fz);
		if( gsl_blas_dnrm2(Dx)<dx || gsl_blas_dnrm2(fx)<eps ) break;
		}

	gsl_matrix_free(Q);
	gsl_matrix_free(R);
	gsl_vector_free(fx);
	gsl_vector_free(z);
	gsl_vector_free(fz);
	gsl_vector_free(df);
	gsl_vector_free(Dx);
}









/*
def newton_jacobian(f:"function", jacobian:"jacobian", x_start:vector, eps:float=1e-3):
	x=x_start.copy()
	while True :
		fx=f(x)
		J=jacobian(x)
		givens.qr(J)
		Dx = givens.solve(J,-fx)
		s=2
		while True :
			s/=2
			y=x+Dx*s
			fy=f(y)
			if fy.norm()<(1-s/2)*fx.norm() or s<0.02 : break
		x=y; fx=fy;
		if fx.norm()<eps : break
	return x;
*/
