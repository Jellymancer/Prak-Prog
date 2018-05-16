#include<float.h>
#include<limits.h>
#include<stdio.h>

int equal(double, double, double, double);

int main(){
int i;
printf("\nAssignment 1: Find the maximum integers\n");
i=1;
while(i+1>i){i++;}
printf("\nmy max int for the while-loop = %i\n",i);
for(i=0;i<i+1;i++){i=i;}
printf("my max int for the for-loop = %i\n",i);
i=0;
do{
i++;}
while(i<i+1);
printf("my max int for the do while loop = %i\n",i);

int j;
while(j>j-1){j--;}
printf("\nThe minimum integer is = %i\n",j);

float f;
double d;
long double l;
d=1.0; f=1.0; l=1.0;
while(1+d!=1){d/=2;} d*=2; 
while(1+l!=1){l/=2;} l*=2; 
while(1+f!=1){f/=2;} f*=2;

printf("\nThe value of the machine epsilon for the different variable types is:\n");
printf("float = %g, double =%g, long double =%Lg\n",f,d,l);
f=1.0;d=1.0;l=1.0;
for(f=1; 1+f!=1; f/=2){}f*=2;
for(d=1; 1+l!=1; l/=2){}l*=2;
for(l=1; 1+d!=1; d/=2){}d*=2;
printf("Using for loops yields:\n");
printf("float = %g, double =%g, long double =%Lg\n",f,d,l);


f=1.0;d=1.0;l=1.0;
do d/=2; while(1+d!=1);d*=2;
do f/=2; while(1+f!=1);f*=2;
do l/=2; while(1+l!=1);l*=2;
printf("Using do while loops yields:\n");
printf("float = %g, double =%g, long double =%Lg\n",f,d,l);

printf("The value of the machine epsilon from the float.h header for the different variables is:\n");
printf("float =%g, double = %g, long double = %Lg",FLT_EPSILON,DBL_EPSILON,LDBL_EPSILON);



float upsum=0;
float downsum=0;
printf("\n\nAssignment 2:  Harmonic series with float and double\n");
int intmax=INT_MAX/2;
for(int i=1; i<=intmax; i++){
		upsum += 1.0/i;
		downsum += 1.0/(intmax-i+1);}
printf("Sum_up_float = %f\n",upsum);
printf("Sum_down_float is %f\n",downsum);
printf("The difference stems from the finite precision of the float variable\n");
printf("This series is a harmonic series and it diverges\n");

double dupsum=0, ddownsum=0;
for(int i=1; i<=intmax; i++){
		dupsum += 1.0/i;
		ddownsum += 1.0/(intmax-i+1);}
printf("Sum_up_double = %g\n",dupsum);
printf("Sum_down_double is %g\n",ddownsum);
printf("This result is much closer since the 'double' variable has twice as many digits\n"); 


printf("\n\nAssignment 3: relative precision\n");
double a=3, b=3.004, tau=0.1, epsilon=0.01;
int ret=equal(a,b,tau,epsilon);
printf("Using a=%g, b=%g, tau=%g, epsilon=%g, the precision function returns: %i\n",a,b,tau,epsilon,ret);
a=5; b=10;
ret=equal(a,b,tau,epsilon);
printf("Using a=%g, b=%g, tau=%g, epsilon=%g, the precision function returns: %i\n",a,b,tau,epsilon,ret);
return 0;}
