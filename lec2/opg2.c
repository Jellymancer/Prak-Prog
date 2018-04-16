#include<stdio.h>


int main(){

float fx = 0.1111111111111111111111111111;
double dx = 0.1111111111111111111111111111;
long double ldx = 0.1111111111111111111111111111L;

printf("\n Assignment 2: finding maximal digits\n");
printf("Printing a long number x=0.1111... using: double, float and long double gives:\n");
printf("float: x=%.25g\n",fx);
printf("double: x=%.25lg\n",dx);
printf("long double: x=%.25Lg\n",ldx);

printf("The place where the digits differ from 1 marks the maximum size of the variable\n");
return 0;
}
