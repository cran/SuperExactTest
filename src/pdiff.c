#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mvhyper.h"
//pdf of set difference (A^B)\C
void C_pdf_ABdiffC(int *x, int *L, int *n, double *p, int *logp){
/*
x:     size of set difference (A^B)\C
L:     subset sizes
n:     background size
p:     output probability
logp:  return log probability
*/
	const double tiny = 1.0E-320;
	int j;
	int i0=0;
	double temp1,temp2;
	int maxJ=min2(L[0],L[1]) - *x;
	int minJ=max2(L[0]+L[1] - *x -*n,0);
	maxJ=min2(L[2],maxJ);
	//sanity check of input must have been performed in the R warpper function
	j = minJ;
	temp1=C_dhyper(*x + j,L[0],*n - L[0],L[1],i0);
	temp2=C_dhyper(j,*x + j,*n - *x - j, L[2],i0);
	*p = temp1*temp2;
	if( fabs(*p) < tiny ) return;
	for(j = minJ+1; j <= maxJ; j++){
		temp1 = temp1 * ((L[0] - *x - j + 1.0)/(*x + j)) * ((L[1] - *x - j + 1.0) /(*n - L[0] - L[1] + *x + j));
		temp2 = temp2 *  ((*x + j + 0.0)/j) * ((L[2] - j + 1.0)/(*n - *x - j + 1.0));
		*p += temp1 * temp2;
	}
	if (*p > 1) *p = 1.0;
	if ( *p < 0 ) *p = db_xmin;
	if(*logp>0) *p=log(*p);
	return;
}
void C_pdf_ABdiffC_logVal(int *x, int *L, int *n, double *p, int *logp, double *logVal){
/*
x:     size of set difference (A^B)\C
L:     subset sizes
n:     background size
p:     output probability
logp:  return log probability
*/
	const double tiny = 1.0E-320;
	int j;
	int i0=0;
	double temp1,temp2;
	int maxJ=min2(L[0],L[1]) - *x;
	int minJ=max2(L[0]+L[1] - *x -*n,0);
	maxJ=min2(L[2],maxJ);
	//sanity check of input must have been performed in the R wrapper function
	j=minJ;
	temp1=C_dhyper_logVal(*x + j,L[1],*n - L[1],L[0],i0,logVal);
	temp2=C_dhyper_logVal(j,*x + j,*n - *x - j, L[2],i0,logVal);
	*p = temp1*temp2;
	if( fabs(*p) <= tiny ) return;
	for(j=minJ+1; j <= maxJ; j++){
		temp1 = temp1 * ((L[1] - *x - j + 1.0)/(*x + j)) * ((L[0] - *x - j + 1.0) /(*n - L[0] - L[1] + *x + j));
		temp2 = temp2 *  ((*x + j + 0.0)/j) * ((L[2] - j + 1.0)/(*n - *x - j +1));
		*p += temp1 * temp2;
	}
	if (*p > 1) *p = 1.0;
	if ( *p < 0 ) *p = db_xmin;
	if(*logp>0) *p=log(*p);
	return;
}
void C_ABdiffC(int *x, int *L, int *n, double *p, int *lower, int *logp){
/*
x:     size of set difference (A^B)\C
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
	double logVal[*n];
	double Xmean;
	int maxX=min2(L[0],L[1]);
	double pp[maxX];
	for(i=1; i<= *n ; i++){
		logVal[i-1]=log((double)i);
	}
	if(*x == 0){
		C_pdf_ABdiffC_logVal(x, L, n, p, &i0, logVal);
		if(*lower == 0) *p = 1.0 - *p;
		if(*logp>0) *p=log(*p);
		return;
	}
	Xmean = 0.0 + *n;
	Xmean = Xmean * ((double) L[0] / *n) * ((double) L[1] / *n) * (1 - (double) L[2] / *n);
	for(i=0; i< maxX ; i++) pp[i]=0.0;
	*p = 0.0;
	if((double) *x > Xmean){
		i = *x + 1;
		for(; i <= maxX; i++){
			C_pdf_ABdiffC_logVal(&i, L, n, &p0, &i0, logVal);
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
			C_pdf_ABdiffC_logVal(&i, L, n, &p0, &i0, logVal);
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
