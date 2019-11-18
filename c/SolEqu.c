#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define RAND_MAX 2147483647
#define PI 3.1415926
#define EPSILON 0.0000000001
#define EPSILON2 0.00001


int ParaLU(double mean, double variance, double * mu, double * sigma);

int main(void){
	
	double mean=3000, variance=1000000, tau, sigma;
	ParaLU(mean, variance, &tau, &sigma);
	printf("tau is %lf sigma is %lf ", tau, sigma);
	
	return 0;}





int ParaLU(double mean, double variance, double * mu, double * sigma){
  double lambda, delta, Max, Min, max, min, xn=2, fx, dfx, error=pow(10,-10);
   
    fx= (xn+1)/(xn-1)*log(xn)- 2*(variance+ mean*mean)/ mean/mean;
    
   while( fabs(fx) >error ) {	  
	  dfx=-2/(xn-1)/(xn-1)* log(xn)+ ( 1+2/(xn-1) )/xn;
	  xn=xn-fx/dfx;	
	  
	  if(xn<=1)
	  xn=(1.0+rand())/(1.0+RAND_MAX)+ 1;
	  
	  fx= (xn+1)/(xn-1)*log(xn)- 2*(variance+ mean*mean)/ mean/mean;	    
	  }
	  
	  lambda=xn;   //lambda=Max/Min, delta= Max-Min
	  delta=mean*log(lambda);
	  Max=lambda/(lambda-1)*delta;
	  Min=delta/(lambda-1);
	  max=log(Max);
	  min=log(Min);
	  *mu= (max+min)/2;
	  *sigma=(max-min)/2;
	  
	  printf("fx %lf, lambda %lf delta %lf Max %lf Min %lf max %lf min %lf mu is %lf sigma is %lf\n", fx, lambda, delta, Max, Min, max, min, *mu, *sigma);
	  
	  return 0;
   }
