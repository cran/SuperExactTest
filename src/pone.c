#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mvhyper.h"
//pdf of set difference A \ union(B,C)
void C_pdf_AdiffBC(int *x, int *L, int *n, double *p, int *logp){
/*
x:     size of set difference A \ union(B,C)
L:     subset sizes
n:     background size
p:     output probability
logp:  return log probability
*/
//	const double tiny = 1.0E-320;
	int jk,k;
	int i0=0;
	double temp1,temp2,tao;
	int maxJK,minJK,maxK,minK;
	minJK=max2(0,L[0] + L[1] - *n);
	maxJK=min2(L[0] - *x,L[1]);
	//sanity check of input must have been performed in the R warpper function
	*p=0.0;
	tao=C_logChoose(*n, L[2]);
	temp1=temp2=0.0;
	for(jk = minJK; jk <= maxJK; jk++){
		if(jk==minJK){
			temp1 = C_dhyper(jk,L[0],*n - L[0],L[1],i0);
		}else{
			temp1 = temp1 * ((L[0]-jk+1.0)/(jk+0.0)) * ((L[1]-jk+1.0) /(*n -L[0]-L[1]+jk));
		}
		maxK = max2(0,min2(L[2] - L[0] + jk + *x,jk));
		minK = max2(L[2] + jk + *x - *n,0);
		for(k = minK; k <= maxK ; k++){
			if(k == minK){
				temp2 = exp(C_logChoose(jk, k)+C_logChoose(L[0]-jk, *x)+C_logChoose(*n - L[0], L[2] - L[0] + jk - k + *x)-tao);
			}else{
				temp2 = temp2 * (jk-k+1.0)/(k+0.0) * (L[2] - L[0] + jk + *x - k + 1.0) / (*n - L[2] - jk + k - *x);
			}
			//Rprintf("%d %d %.6e %.6e\n",jk,k,temp1,temp2);
			*p += temp1 * temp2;
		}
	}
	if (*p > 1) *p = 1.0;
	if ( *p < 0 ) *p = db_xmin;
	if(*logp>0) *p=log(*p);
	return;
}
void C_AdiffBC(int *x, int *L, int *n, double *p, int *lower, int *logp){
/*
x:     size of set difference A \ union(B,C)
L:     subset sizes
n:     background size
p:     output probability
lower: 1, lower tail probability Pr(overlap <= x); 0, upper tail probability Pr(overlap > x)
logp:  return log probability
*/
	const double tiny = 1.0E-320;
	int i,j;
	int i0=0;
	double p0=0.0;
	double Xmean;
	int maxX;
	maxX=L[0];
	double pp[maxX];
	if(*x == 0){
		C_pdf_AdiffBC(x, L, n, p, &i0);
		if(*lower == 0) *p = 1.0 - *p;
		if(*logp>0) *p=log(*p);
		return;
	}
	Xmean = L[0] * (1.0 + ((L[1]+0.0) / *n) * ((L[2]+0.0) / *n) - ((L[1]+L[2]+0.0) / *n));
	for(i=0; i< maxX ; i++) pp[i]=0.0;
	*p = 0.0;
	if((double) *x > Xmean){
		i = *x + 1;
		for(; i <= maxX; i++){
			C_pdf_AdiffBC(&i, L, n, &p0, &i0);
			pp[i]=p0;
			if(p0 <= tiny) break;
			if(i > (*x + 1) && (p0/pp[i-1]) <= 0.01) break;  //stop if no improve in precision
		}
		for(j = i; j >= *x + 1; j--){ //iteration from smallest to largest; more accurate
			*p += pp[j];
		}
		if(*lower > 0) *p = 1.0 - *p;
	}else{
		i = *x;
		for(;i >= 0; i--){
			C_pdf_AdiffBC(&i, L, n, &p0, &i0);
			pp[i]=p0;
			if(p0 <= tiny) break;
			if(i < *x && (p0/pp[i+1]) <= 0.01) break;  //stop if no improve in precision
		}
		for(j = i; j <= *x; j++){ //iteration from smallest to largest; more accurate
			*p += pp[j];
		}
		if(*lower == 0) *p = 1.0 - *p;
	}
	if (*p > 1) *p = 1.0;
	if ( *p < 0 ) *p = db_xmin;
	if(*logp > 0) *p = log(*p);
	return;
}
