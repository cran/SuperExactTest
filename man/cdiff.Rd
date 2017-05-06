\name{dpdiff}
\alias{dpdiff}
\alias{cpdiff}
\alias{dpone}
\alias{cpone}
\docType{package}
\title{Set Difference Probability}
\description{
  Density and distribution functions of set difference tests.
}
\usage{
dpdiff(x,s1,s2,s3,n,log.p=FALSE)
dpone(x,s1,s2,s3,n,log.p=FALSE)
cpdiff(x,s1,s2,s3,n,lower.tail=TRUE,log.p=FALSE)
cpone(x,s1,s2,s3,n,lower.tail=TRUE,log.p=FALSE)
}
\arguments{
 \item{x}{integer, number of elements that are different.}
 \item{s1,s2,s3}{set sizes.} 
 \item{n}{integer, background population size.}
 \item{lower.tail}{logical; if TRUE, probability is \code{P[X <= x]}, otherwise, \code{P[X > x]}.}
 \item{log.p}{logical; if \code{TRUE}, probability p is given as \code{log(p)}.}
}
\value{
\code{dpdiff} gives the density and \code{cpdiff} gives the distribution function for set difference of \code{intersect(A,B) \\ C}, while \code{dpone} gives the density and \code{cpone} gives the distribution function of set difference of \code{A \\ union(B,C)}, where \code{|A|=s1}, \code{|B|=s2}, and \code{|C|=s3}.
}
\author{
Minghui Wang <minghui.wang@mssm.edu>, Bin Zhang <bin.zhang@mssm.edu>
}
\references{
Minghui Wang, Yongzhong Zhao, and Bin Zhang (2015). Efficient Test and Visualization of Multi-Set Intersections. \emph{Scientific Reports} 5: 16923.
}
\seealso{\code{\link{cpsets}}, \code{\link{dpsets}}}
\examples{
\dontrun{
#set up fake data
n=500; s1=260; s2=320; s3=430; x=80
(d1=sapply(0:80,function(x) SuperExactTest:::dpdiff(x,s1,s2,s3,n)))
(p1=sapply(0:80,function(x) SuperExactTest:::cpdiff(x,s1,s2,s3,n,lower.tail=TRUE)))
par(mfrow=c(1,2))
plot(0:80,d1,main='Density');
plot(0:80,p1,main='Distribution');
}
}
