#include<stdio.h>
#include<math.h>

double x;
int main()
{
while( scanf("%lg",&x) != EOF) {printf("%lg \t %lg\n",x,cos(x));}
return 0;
}
