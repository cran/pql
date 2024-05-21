#' The average weighted estimator and the unknown weighted estimator of the PQL in #'Poisson-GLMS through damped Gauss-Newton
#' @param data is a design matrix with uniform distribution and the response vector
#' @param G is the number of subsets.
#' nk is the size of subsets.
#' n is the sample size.
#' p is the number of variables.
#' @return betaBA, betaBW, MSEA, MSEW
#' @export
#' @examples
#'p<- 5;G<- 20;n<- 1000;nk=200
#'X<- matrix(runif(n * p, 0, 0.5), ncol = p)
#'beta <- runif(p, 0, 1)
#'y<- rpois(n, exp(X%*% beta))
#'data=cbind(y,X)
#'  pqlBpoisson2(data,G,nk)
pqlBpoisson2=function(data,G,nk){
  betaA=pqlPoisson(data,G,nk)$betaA
  gbetap= y%*%X-t(t(X)%*%exp(X%*%betaA))
  gbetap#SQ
  p=5
  n= 1000
  X<- matrix(runif(n * p, 0, 0.5), ncol = p)
  beta <- runif(p, 0, 1)
  y<- rpois(n, exp(X %*% beta))
  hbetap=matrix(rep(0,p*p),ncol=p)
  hbetap= t(X) %*%X
  hbetap;hA=hbetap;gA=gbetap
  betaBA=betaA+t(solve(hA)%*%t(gA))/G
  betaW=pqlPoisson(data,G,nk)$betaW
  gbetap= y%*%X- t(t(X)%*%exp(X%*%betaW))
  hbetap= t(X) %*%X
  hbetap;hw=hbetap;gw=gbetap
  betaBW=betaW-t(solve(hw)%*%t(gw))/G
  beta <- as.numeric(beta)
  betaBA <- as.numeric(betaBA)
  norm(as.matrix(beta-betaBA),"f");norm(as.matrix(beta-betaBW),"f")
  return(list(betaBA=betaBA,betaBW=betaBW,MSEA=norm(as.matrix(beta-betaBA),"f"),
              MSEW=norm(as.matrix(beta-betaBW),"f")))
}

