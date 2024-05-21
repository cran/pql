#' A partitioned quasi-likelihood for distributed statistical inference
#' @param data is a design matrix with uniform distribution and the response vector
#' @param G is the number of subsets.
#'        nk is the size of subsets.
#'        n is the sample size.
#'        p the number of variables.
#' @return betaBA, betaBW, MSEA, MSEW
#' @export
#' @examples
#'library(parallel)
#'library(numDeriv)
#'library(Rmpi)
#'install.packages("pracma");
#'library(pracma)
#'p<- 5;G<- 20;n<- 1000;nk=200
#'X<- matrix(runif(n * p, 0, 0.5), ncol = p)
#'beta <- runif(p, 0, 1)
#'y<- rpois(n, exp(X%*% beta))
#'data=cbind(y,X)
#'  pqlBpoisson1(data,G,nk)
pqlBpoisson1=function(data,G,nk){
  beta <- runif(p, 0, 1)
  p=ncol(data[,-1]);n=nrow(data[,-1])
  y=data[,1];X=data[,-1]
  betaA=pqlPoisson(data,G,nk)$betaA
  l= function(beta){
    lbeta=y%*% X%*% beta-sum(exp(X %*% beta))-sum(log(prod(y[!y==0])))
    lbeta}
  gA=grad(l, betaA)
  hA=hessian(l, betaA)
  betaBA=betaA+t(solve(t(hA)%*%hA,t(hA)%*%gA))/G
  betaW=pqlPoisson(data,G,nk)$betaW
  gW=grad(l, betaW)
  hW=hessian(l, betaW)
  beta <- as.numeric(beta)
  betaBA <- as.numeric(betaBA)
  betaBW=betaW-t(solve(t(hW)%*%hW,t(hW)%*%gW)) /G
  norm(as.matrix(beta-betaBA),"f");norm(as.matrix(beta-betaBW),"f")
  return(list(betaBA=betaBA,betaBW=betaBW,MSEA=norm(as.matrix(beta-betaBA),"f"),
              MSEW=norm(as.matrix(beta-betaBW),"f")))
}
