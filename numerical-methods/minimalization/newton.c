#include "func.h"

int newton(double f(gsl_vector* x,gsl_vector* gr_fx), void Hess(gsl_vector* p, gsl_matrix* H),gsl_vector* x, double dx, double eps){
	int n=x->size;
	gsl_matrix* H = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* gr_fx = gsl_vector_alloc(n);
	gsl_vector* z  = gsl_vector_alloc(n); // z used for bactracking linesearch.
	gsl_vector* gr_fz = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);
	int ncalls = 0;
	do{ ncalls ++; f(x,gr_fx);

		Hess(x,H); //Set Hessian

		qrgsdecomp(H,R);
		qrgssolve(H,R,gr_fx,Dx);
		gsl_vector_scale(Dx,-1); //finding dx step by solving df(x)+H(x)*dx by QR decomposition.
		double lam=1,fvalz, fvalx, armcond; double *armcondp=&armcond;

		do{
			fvalz=f(z,gr_fz);
			fvalx=f(x,gr_fx);
			gsl_blas_ddot(Dx,gr_fx,armcondp); armcond = *(armcondp);

			lam/=2.0;
			gsl_vector_memcpy(z,x); //Bactracking linesearch is performed until condition armijo condition is met.
			gsl_blas_daxpy(0.5,Dx,z);

			} while(fvalz > fvalx + 1e-4*lam*armcond && lam >0.02);
			//Stop when armijo condition is met or lambda becomes too small.

		gsl_vector_memcpy(x,z);
		gsl_vector_memcpy(gr_fx,gr_fz);
		}while(gsl_blas_dnrm2(Dx)>dx && gsl_blas_dnrm2(gr_fx)>eps);

	return ncalls;
	gsl_matrix_free(H);
	gsl_matrix_free(R);
	gsl_vector_free(gr_fx);
	gsl_vector_free(z);
	gsl_vector_free(gr_fz);
	gsl_vector_free(Dx);
}

int newton_broy(double f(gsl_vector* x,gsl_vector* gr_fx),gsl_vector* x, double eps){
	int n=x->size;
	gsl_matrix* H = gsl_matrix_alloc(n,n); // Here H is the inverse hessian. Updates are performed directly on H⁻¹.
	gsl_vector* bl1 = gsl_vector_alloc(n); // for blas routines
	gsl_vector* bl2 = gsl_vector_alloc(n); //
	gsl_vector* gr_fx = gsl_vector_alloc(n);
	gsl_vector* z  = gsl_vector_alloc(n); // z used for bactracking linesearch.
	gsl_vector* gr_fz = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n); //step

	int ncalls = 0, callmax=1e8;
	gsl_matrix_set_identity(H); //Set hessian as identity (first approx)
	f(x,gr_fx);
	double lam,fvalz,fvalx,armcond,Hys, alpha = 1e-4;
	double *armcondp=&armcond; double *Hysp=&Hys;

		do{ncalls++;
		lam = 1;
		//Solve  grad_f(x + s) =  grad_f(x) + (H+dH)*s, by setting the gradient to zero.
		// s=-H⁻¹ grad_f(x)
		f(x,gr_fx);
		gsl_blas_dgemv(CblasNoTrans,-1.0,H,gr_fx,0.0,Dx); // set the step.

		while(1){
                        gsl_vector_memcpy(z,x); //Bactracking linesearch is performed until condition armijo condition is met.
                        gsl_vector_add(z,Dx);//adds the step to x.

                        fvalz=f(z,gr_fz);
                        fvalx=f(x,gr_fx);
                        gsl_blas_ddot(Dx,gr_fx,armcondp); armcond = *(armcondp);
			if(lam < 0.002){ gsl_matrix_set_identity(H); break;} //if update diverges, reset Hessian.
			if(fvalz < fvalx + alpha*lam*armcond) break;
			lam/=2.0;
			gsl_vector_scale(Dx,0.5);
                        }                        //Stop when armijo condition is met or lambda becomes too small.

		f(z,gr_fz);
		gsl_vector_sub(gr_fz,gr_fx);
		// Make the broyden update. H⁻¹ -> H⁻¹ + (s-H⁻¹y)*s^T*H⁻¹ / (y^T*H⁻¹ s)
		gsl_blas_dgemv(CblasNoTrans,1.0,H,gr_fz,0.0,bl1); // H⁻¹y ->bl1
		gsl_blas_dgemv(CblasNoTrans,1.0,H,Dx,0.0,bl2); //H⁻¹s ->bl2
		gsl_blas_ddot(bl1,Dx,Hysp); Hys=*(Hysp); Hys=1/Hys; // 1/yH⁻¹s->Hys
		gsl_vector_sub(Dx,bl1); // s-H⁻¹y->s
		gsl_blas_dger(Hys,Dx,bl2,H);//updating the H⁻¹

		gsl_vector_memcpy(x,z);
		f(x,gr_fx);
		if(ncalls >= callmax) printf("ERROR:MAX ITERATIONS ACHIEVED!");
		}while(gsl_blas_dnrm2(gr_fx)>eps && ncalls < callmax);

	return ncalls;
	gsl_matrix_free(H);
	gsl_vector_free(bl1);
	gsl_vector_free(bl2);
	gsl_vector_free(gr_fx);
	gsl_vector_free(z);
	gsl_vector_free(gr_fz);
	gsl_vector_free(Dx);
}

