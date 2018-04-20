#include<stdio.h>
#include "nvector.h"


nvector* nvector_alloc(int n){
  nvector* v = malloc(sizeof(nvector));
  (*v).size = n;
  (*v).data = malloc(n*sizeof(double));
  if( v==NULL ) fprintf(stderr,"error in allocation\n");
  return v;
}

void nvector_free(nvector* v){ free(v->data); free(v);}

void nvector_set(nvector* v, int i, double value){ (*v).data[i]=value; }

double nvector_get(nvector* v, int i){return (*v).data[i]; }


void nvector_print(char* s, nvector* v){
	printf("%s",s);
	for(int i =0;i < v->size;i++){
		printf("%g\n",(*v).data[i]);
	}}


void nvector_set_zero(nvector* v){
	for(int i=0;i<v->size;i++) nvector_set(v,i,0);}


int double_equal(double a, double b){
	double TAU = 1e-6, EPS = 1e-6;
	if (fabs(a - b < TAU))
		return 1;
	if (fabs(a - b) / (fabs(a) + fabs(b)) < EPS / 2)
		return 1;
	return 0;}

int nvector_equal(nvector* a, nvector* b){
	if (a->size != b->size) return 0;
	for (int i = 0; i < b->size; i++)
		if (!double_equal(b->data[i], a->data[i]))
			return 0;
	return 1;}

void nvector_add(nvector* a, nvector* b,nvector* c){
	if(a->size != b->size) fprintf(stderr,"Sizes of vectors are not equal\n");
	if(a->size != c->size) fprintf(stderr,"Sizes of vectors are not equal\n");
	if(b->size != c->size) fprintf(stderr,"Sizes of vectors are not equal\n");

	for(int i=0; i<a->size;i++){
		nvector_set(c,i,(*a).data[i]+(*b).data[i]);}}

void nvector_sub(nvector* a, nvector* b,nvector* c){
	if(a->size != b->size) fprintf(stderr,"Sizes of vectors are not equal\n");
	if(a->size != c->size) fprintf(stderr,"Sizes of vectors are not equal\n");
	if(b->size != c->size) fprintf(stderr,"Sizes of vectors are not equal\n");

	for(int i=0; i<a->size;i++){
		nvector_set(c,i,(*a).data[i]-(*b).data[i]);}}


void nvector_scale(nvector* v, nvector* c, double x){
	if(v->size != c->size) fprintf(stderr,"Sizes of vectors are not equal\n");
	for(int i=0;i<v->size;i++) nvector_set(c,i,nvector_get(v,i)*x);}
