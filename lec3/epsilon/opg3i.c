#include<stdio.h>


int equal(double a,double b,double tau,double epsilon)
{
int r;

if(a>=b && (a-b)<tau){return 1;}
else if(a<b && (b-a)<tau){return 1;}
else if(a<b && ((b-a) / (a+b)) <= (epsilon/2)){return 1;}
else if(a>b && ((a-b) / (a+b)) < (epsilon/2)){return 1;}
else{return 0;}

}


