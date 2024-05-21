\name{pqlBpoisson2}
\alias{pqlBpoisson2}
\title{The weighted Gauss-Newton estimators of the PQL in Poisson-GLMS}
\usage{
pqlBpoisson2(data,G,nk)}
\arguments{
\item{data}{is a design matrix with uniform distribution and the response vector}
\item{G}{is the number of subsets.}
\item{nk}{is the size of subsets.}
}
\value{
betaBA, betaBW, MSEA, MSEW}
\description{
The average weighted estimator and the unknown weighted estimator of the PQL in Poisson-GLMS through damped Gauss-Newton}
\examples{
p<- 5;G<- 20;n<- 1000;nk=50
X<- matrix(runif(n * p, 0, 0.5), ncol = p)
beta =matrix(runif(p, 0, 1),nrow=p)
L=X\%*\%beta
y<- rpois(n, exp(L))
data=cbind(y,X)
pqlBpoisson2(data,G,nk)
}

