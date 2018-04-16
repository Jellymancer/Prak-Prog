#include <stdio.h>
#include <complex.h>
#include <math.h>
int main()
{
double valuer,valuei;
double complex result;

valuer=0;
valuei=1;
result = cexp(valuer+valuei*I);
printf("\nNow I calcualte the complex valued problems:\n");
printf("Complex exponential of  %f%+fi is %f%+fi\n",valuer,valuei,creal(result),cimag(result));

valuer=0;
valuei=1;
result=cexp(valuer*M_PI+valuei*I*M_PI);
printf("Complex exponmential of %f%+fi is %.2f+%.2fi\n  ",valuer*M_PI,valuei*M_PI,creal(result),cimag(result));

valuer=0;
valuei=1;
result=cpow(valuer+valuei*I,M_E);
printf("%f%+fi to the power of e is %f%+.2fi\n",valuer,valuei,creal(result),cimag(result));
return 0;
}
