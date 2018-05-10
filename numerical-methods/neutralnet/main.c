#include "neurons.h"

double activation_func(double x){return x*exp(-x*x);}
double fit_func(double x){return cos(5*x-1)*exp(-x*x);} // functions for fig 1

double activation_func2(double x){return cos(x)*exp(-x*x);} // and fig 2.
double fit_func2(double x){return x*x*x;}



int main(){
//Fig 1. Interpolating wavelet using x*exp(-x*x)
	int n=5; // two neurons
	neurons* nw=neurons_alloc(n,activation_func);
	double a=-1,b=1; //x-limits.
	int nx=20; // allocate list of x and y for interpolation and fill with values.
	gsl_vector* vx=gsl_vector_alloc(nx); //x-list
	gsl_vector* vy=gsl_vector_alloc(nx); // y-list
	for(int i=0; i<nx;i++){
		double x=a+(b-a)*i/(nx-1);
		double f = fit_func(x);
		gsl_vector_set(vx,i,x);
		gsl_vector_set(vy,i,f);
}

//setting starting parameters ai,bi and weigths w.
	for(int i=0;i<nw->n;i++){
	gsl_vector_set(nw->data,0*nw->n+i,a+(b-a)*i/(nw->n-1)); //ai
	gsl_vector_set(nw->data,1*nw->n+i,1); //bi
	gsl_vector_set(nw->data,2*nw->n+i,1); //wi
}
	neurons_train(nw,vx,vy); //train the network
//printing the answers
	FILE* nn = fopen("nn_out.txt","w");
	for(int i=0;i<vx->size;i++){
		fprintf(nn,"%g %g\n",gsl_vector_get(vx,i),gsl_vector_get(vy,i));
	}
	fprintf(nn,"\n\n");
	double dz=0.001;
	for(double z=a; z<=b;z+=dz){
		double y=neurons_feed_forward(nw,z);
		fprintf(nn,"%g %g\n",z,y);
	}
	fclose(nn);



// Fig 2. Interpolating x³ using wavelet.
	n=5; // no of neurons
	neurons* nw2=neurons_alloc(n,activation_func2);
	a=-2, b=2;
	nx=20; // allocate list of x and y for interpolation and fill with values.
	gsl_vector* vx2=gsl_vector_alloc(nx); //x-list
	gsl_vector* vy2=gsl_vector_alloc(nx); // y-list

	for(int i=0; i<nx;i++){
		double x=a+(b-a)*i/(nx-1);
		double f = fit_func2(x);
		gsl_vector_set(vx2,i,x);
		gsl_vector_set(vy2,i,f);
}

	for(int i=0;i<nw2->n;i++){
	gsl_vector_set(nw2->data,0*nw2->n+i,a+(b-a)*i/(nw2->n-1)); //ai
	gsl_vector_set(nw2->data,1*nw2->n+i,1); //bi
	gsl_vector_set(nw2->data,2*nw2->n+i,1); //wi
}
 	neurons_train(nw2,vx2,vy2);
	FILE* nn2 = fopen("nn_out2.txt","w");
	for(int i=0;i< vx2->size;i++){
		fprintf(nn2,"%g %g\n",gsl_vector_get(vx2,i),gsl_vector_get(vy2,i));
	}
	fprintf(nn2,"\n\n");
	for(double z=a; z<=b;z+=dz){
		double y=neurons_feed_forward(nw2,z);
		fprintf(nn2,"%g %g\n",z,y);
	}
	fclose(nn2);

printf("\nANN interpolation is performed for two types of tabluated using two different activation functions.\n");
printf("Fig 1: shows the interpolation of f(x)=cos(5*x-1)*exp(-x²) using activation function x*exp(-x²)\n");
printf("Fig 2: shows the interpolation of f(x)=x³ using activation function cos(x)*exp(-x²)\n");




neurons_free(nw);
neurons_free(nw2);
return 0;
}
