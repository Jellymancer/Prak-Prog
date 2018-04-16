#include <stdio.h>
#include <math.h>
int main()
{
float result,value;
value=5;
result = tgamma(value);
printf("\nAssignment 1: Mathematical functions.\n First the real values functions:\n");
printf("Gamma function of %.2f is %.2f  \n",value,result);

value=0.5;
result=j0(value);
printf("Bessel function of %.2f is %.2f  \n",value,result);

value=-2;
result=sqrt(value);
printf("Square root of %.2f is %.2f \n",value,result);
return 0;
}
