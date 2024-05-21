\name{ pqlBpoisson1}
\alias{ pqlBpoisson1}
\title{The weight Gauss-Newton estimators of the PQL in Poisson-GLMS}
\usage{
pqlBpoisson1(data,G,nk)}
\arguments{
\item{data}{is a design matrix with uniform distribution and the response vector}
\item{G}{is the number of subsets.}
\item{nk}{is the size of subsets}

}
\value{betaBA, betaBW, MSEA, MSEW}
\description{The average weighted estimator and the unknown weighted estimator of the PQL in Poisson-GLMS through damped Gauss-Newton
}
\examples{
G <- 20;n=1000;p=5; nk=50
X<- matrix(runif(1000* 5, 0, 0.5), ncol = 5)
beta =matrix(runif(p, 0, 1),nrow=p)
L=X\%*\%beta
y<- rpois(1000, exp(L))
data=cbind(y,X)
pqlBpoisson1(data,G,nk)}

