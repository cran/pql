#' The average weighted estimator and the unknown weighted estimator of the PQL in #' Poisson-GLMS through damped Gauss-Newton
#' @description A partitioned quasi-likelihood for distributed statistical inference
#' @param data is a design matrix with uniform distribution and the response vector
#' @param G is the number of subsets.
#' nk is the number of outer subsets.
#' n is the sample size.
#' p is the number of variables.
#' @return betaBA, betaBW, MSEA, MSEW
#' @export
#' @examples
#'library(parallel)
#'library(numDeriv)
#'library(Rmpi)
#'install.packages("pracma");
#'library(pracma)
#' p<- 5;G<- 20;n<- 1000;nk=200
#'X<- matrix(runif(n * p, 0, 0.5), ncol = p)
#' beta <- runif(p, 0, 1)
#'y<- rpois(n, exp(X%*% beta))
#'data=cbind(y,X)
#' pqlPoisson(data,G,nk)

pqlPoisson=function(data,G,nk){
  p=ncol(data[,-1]);n=nrow(data[,-1])
  y=data[,1];X=data[,-1]
  b <- runif(p, 0, 1)
  beta =matrix(b,nrow=p)
  be=matrix(rep(0,G*p),ncol=p);v=w=rep(0,G)
  Rm=matrix(rep(0, nk*G),ncol=G)
  mr=matrix(rep(0,nk*G),ncol=nk)
  R=matrix(rep(0,nk*n),ncol=n)
  for(i in 1:G )  {
    mr[i,]=sample(1:n,nk,replace=F)
    r=matrix(c(1:nk,mr[i,],ncol=nk,byrow=T))
    R[t(r)]=1
    Xr= R%*%X
    yr=R%*%y
    be[i,]<-glm(yr~ Xr-1, family = quasipoisson)$coefficients
    #be[i,]=coef(glm(yr~ Xr-1, family = quasipoisson)
    v[i]=1/var(be[i,])
    w[i]=(v[i])/(sum(v))
  }
  betaA=rep(0,p);betaW=rep(0,p)
  for (j in 1:p) {
    betaA[j]=sum(be[,j])/G
    betaW[j]=t(be[,j])%*%w
  }
  beta <- as.numeric(beta)
  betaA <- as.numeric(betaA)
  betaW <- as.numeric(betaW)
  var(beta-betaA);var(beta-betaW)
  return(list(betaW=betaW,betaA=betaA,MSEW=var(beta-betaW),MSEA=var(beta-betaA)))
}

