#Compute multi set intersection test
#Author: Minghui Wang, minghui.wang@mssm.edu
#Date: 20 July, 2014

#Input:
#x           number of overlap elements or set difference between the subsets
#L           a vector of subset sizes
#n           background size
#log.p       logical; if TRUE, probability p is given as log(p).
#lower.tail  logical; if TRUE (default), probability is P[overlap <= x], otherwise, P[overlap > x].

#distribution function
cpsets <- function(x,L,n,lower.tail=TRUE,log.p=FALSE,simulation.p.value=FALSE,number.simulations=1000000){
	if(length(L)<2) stop('L should have at least two entries\n')
	if(x<0 || n<1 || any(L>n) || any(L<x)){
		stop('Invalid input\n')
	}
	if(simulation.p.value) return(cpsets.simulation(x,L,n,lower.tail,log.p,number.simulations))
	L=sort(L)
	.C("C_pmvhyper",as.integer(x),length(L),as.integer(L),as.integer(n),as.numeric(0.0),as.integer(lower.tail),as.integer(log.p))[[5]]
}
cpsets.simulation <- function(x,L,n,lower.tail=TRUE,log.p=FALSE,number.simulations=1000000){
	nL=length(L)
	cc=sapply(1:number.simulations,function(i){
		Ls=sapply(L,function(m) sample.int(n=n,size=m,replace=FALSE))
		sum(table(unlist(Ls))==nL) <= x
	})
	p=sum(cc)/number.simulations
	if(lower.tail==F) p=1-p
	ifelse(log.p,log(p),p)
}
cpdiff <- function(x,s1,s2,s3,n,lower.tail=TRUE,log.p=FALSE){
	L=c(s1,s2,s3)
	if(x<0 || n<1 || any(L>n) || any(L[1:2]<x)) stop('Invalid input\n')
	.C("C_ABdiffC",as.integer(x),as.integer(L),as.integer(n),as.numeric(0.0),as.integer(lower.tail),as.integer(log.p))[[4]]
}
#density function
dpsets <- function(x,L,n,log.p=FALSE){
	if(x<0 || any(L>n) || any(L<x)) if(!log.p){return(0)}else{stop('Invalid input\n')}
	nL=length(L)
	if(nL<2) stop('L should have at least two entries\n')
	L=sort(L)
	.C("C_dmvhyper",as.integer(x),length(L),as.integer(L),as.integer(n),as.numeric(0.0),as.integer(log.p))[[5]]
}
dpdiff <- function(x,s1,s2,s3,n,log.p=FALSE){
#probability density of | intersect(A,B) \ C | = x, where |A|=s1, |B|=s2, |C|=s3
	L=c(s1,s2,s3)
	if(x<0 || any(L>n) || any(L[1:2]<x)) if(!log.p){return(0)}else{stop('Invalid input\n')}
	.C("C_pdf_ABdiffC",as.integer(x),as.integer(L),as.integer(n),as.numeric(0.0),as.integer(log.p))[[4]]
}
##### dpdiff2 and dpdiff3 are R purely versions of dpdiff
##dpdiff2 <- function(x,s1,s2,s3,n,log.p=FALSE){
###probability density of | intersect(A,B) \ C | = x, where |A|=s1, |B|=s2, |C|=s3
##	L=c(s1,s2,s3)
##	if(x<0 || any(L>n) || any(L[1:2]<x)) if(!log.p){return(0)}else{stop('Invalid input\n')}
##	ltao=lchoose(n,s2)+lchoose(n,s3)
##	p=0
##	for(j in seq(0,min(s1-x,s2-x,s3),by=1)){
##		p=p+exp(lchoose(s1,j+x)+lchoose(n-s1,s2-j-x)+lchoose(j+x,j)+lchoose(n-j-x,s3-j)-ltao)
##	}
##	if(log.p) p=log(p)
##	p
##}
##dpdiff3 <- function(x,s1,s2,s3,n,log.p=FALSE){
###probability density of | intersect(A,B) \ C | = x, where |A|=s1, |B|=s2, |C|=s3
##	L=c(s1,s2,s3)
##	if(x<0 || any(L>n) || any(L[1:2]<x)) if(!log.p){return(0)}else{stop('Invalid input\n')}
##	p=0
##	for(j in seq(max(0,s1+s2-x-n),min(s1-x,s2-x,s3),by=1)){
##		p=p+dhyper(j+x,s1,n-s1,s2)*dhyper(j,j+x,n-j-x,s3)
##	}
##	if(log.p) p=log(p)
##	p
##}
cpone <- function(x,s1,s2,s3,n,lower.tail=TRUE,log.p=FALSE){
	L=c(s1,s2,s3)
	if(x<0 || any(L>n) || any(L<0) || x > s1) if(!log.p){return(0)}else{stop('Invalid input\n')}
	.C("C_AdiffBC",as.integer(x),as.integer(L),as.integer(n),as.numeric(0.0),as.integer(lower.tail),as.integer(log.p))[[4]]
}
dpone <- function(x,s1,s2,s3,n,log.p=FALSE){
#probability density of | A \ union(B,C) | = x, where |A|=s1, |B|=s2, |C|=s3
	L=c(s1,s2,s3)
	if(x<0 || any(L>n) || any(L<0) || x > s1) if(!log.p){return(0)}else{stop('Invalid input\n')}
	.C("C_pdf_AdiffBC",as.integer(x),as.integer(L),as.integer(n),as.numeric(0.0),as.integer(log.p))[[4]]
}
### dpone2 and dpone3 are R versions of dpone
###dpone2 <- function(x,s1,s2,s3,n,log.p=FALSE){
####probability density of | A \ union(B,C) | = x, where |A|=s1, |B|=s2, |C|=s3
###	L=c(s1,s2,s3)
###	if(x<0 || any(L>n) || any(L<0) || x > s1) if(!log.p){return(0)}else{stop('Invalid input\n')}
###	ltao=lchoose(n,s2)+lchoose(n,s3)
###	p=0
###	for(j in seq(0,min(s1-x,s2),by=1)){
###		for(k in seq(0,min(s1-x-j,s2-j),by=1)){
###			p=p+exp(lchoose(s1,j+k)+lchoose(n-s1,s2-j-k)+lchoose(j+k,j)+lchoose(s1-j-k,x)+lchoose(n-s1,s3+j+x-s1)-ltao)
###		}
###	}
###	if(log.p) p=log(p)
###	p
###}
###dpone3 <- function(x,s1,s2,s3,n,log.p=FALSE){
####probability density of | A \ union(B,C) | = x, where |A|=s1, |B|=s2, |C|=s3
###	L=c(s1,s2,s3)
###	if(x<0 || any(L>n) || any(L<0) || x > s1) if(!log.p){return(0)}else{stop('Invalid input\n')}
###	tao=lchoose(n,s3)
###	p=0
###	for(jk in seq(0,min(s1-x,s2),by=1)){
###		temp1=dhyper(jk,s1,n-s1,s2)
###		for(k in seq(0,jk,by=1)){
###			temp2=exp(lchoose(jk,k)+lchoose(s1-jk,x)+lchoose(n-s1,s3+jk+x-s1-k)-tao)
###			#cat(jk,k,temp1,temp2,'\n')
###			p=p+temp1*temp2
###		}
###	}
###	if(log.p) p=log(p)
###	p
###}
###cpone2 <- function(x,s1,s2,s3,n,lower.tail=TRUE,log.p=FALSE){
###	p=sapply(seq(0,x,by=1),function(x,s1,s2,s3,n) dpone(x,s1,s2,s3,n),s1=s1,s2=s2,s3=s3,n=n)
###	p=sort(p)
###	p=sum(p)
###	if(!lower.tail) p=1-p
###	if(log.p) p=log(p)
###	p
###}
#
##cpbm <- function(x,A1,A2,B1,B2,n,lower.tail=TRUE,log.p=FALSE){
##	return(phyper(x,A1,n-A1,B1,lower.tail=lower.tail,log.p=log.p))
##	#
##	p=0
##	for(i in seq(0,x,by=1)) p=p+dpbm(i,A1,A2,B1,B2,n)
##	if(!lower.tail) p=1-p
##	if(log.p) p=log(p)
##	p
##}
##dpbm <- function(x,A1,A2,B1,B2,n,log.p=FALSE){
###probability density of the overlap between subsets A1 and B1, A1 and A2 are non-overlappling and B1 and B2 are non-overlapping
###this is essentially the same as dhyper(x,A1,n-A1,B1,log=log.p), i.e., two samples of size A1 and B1 sharing x common elements
##	return(dhyper(x,A1,n-A1,B1,log=log.p))
##	#
##	L=c(A1,B1,A2,B2)
##	if(any(L<0) || any(L>n) || x > A1 || x > B1 || A1+A2 > n || B1+B2>n) if(!log.p){return(0)}else{stop('Invalid input\n')}
##	ltao=lchoose(n,B2)+lchoose(n-B2,B1)
##	p=0
##	for(j in seq(0,min(A2,B1),by=1)){
##		for(k in seq(0,min(A1-x,B2),by=1)){
##			for(l in seq(0,min(A2-j,B2),by=1)){
##				if (n-A1-A2 < B2-k-l) next
##				if (n-A1-A2-B2+k+l < B1-j-x) next
##				p=p+exp(lchoose(A1,k)+lchoose(A1-k,x)+lchoose(A2,l)+lchoose(A2-l,j)+lchoose(n-A1-A2,B2-k-l)+lchoose(n-A1-A2-B2+k+l,B1-j-x)-ltao)
##			}
##		}
##	}
##	if(log.p) p=log(p)
##	p
##}
#Sample usage:
#fake data
#n=500; A=260; B=320; C=430; D=300; x=170; L=c(A,B,C,D)
#(d=dpsets(x,L,n))
#(p=cpsets(x,L,n,lower.tail=F))
