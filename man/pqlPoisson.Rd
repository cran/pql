\name{pqlPoisson }
\alias{pqlPoisson }
\title{The weighted Gauss-Newton estimators of the PQL in Poisson-GLMs}
\usage{
pqlPoisson(data,G,nk)}
\arguments{
\item{data}{is a design matrix with uniform distribution and the response vector}
\item{G}{is the number of subsets}
\item{nk}{is the number of outer subsets.}
}
\value{
betaBA, betaBW, MSEA, MSEW}
\description{
The average weighted estimator and the unknown weighted estimator of the PQL in Poisson-GLMS through damped Gauss-Newton}
\examples{
#library(parallel)
#library(numDeriv)
#library(Rmpi)
#install.packages("pracma");
#library(pracma)
p<- 5;G<- 20;n<- 1000;nk=200
X<- matrix(runif(n*p, 0, 0.5), ncol = p)
beta =matrix(runif(p, 0, 1),nrow=p)
L=X\%*\%beta
y<- rpois(n, exp(L))
data=cbind(y,X)
pqlPoisson(data,G,nk)
}

