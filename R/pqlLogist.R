#' @param data is a total data set
#' @param G is the number of nodes
#' @param nk is  the length of each data subset
#' @return betaW,betaA,MSEW,MSEA
#' @export
#' @examples
#'library(parallel)
#'library(numDeriv)
#'library(Rmpi)
#'install.packages("pracma");
#'library(pracma)
#'G <- 20;n=1000;p=5; nk=50
#'b <- runif(p, 0, 1)
#'beta =matrix(b,nrow=p)
#'X=matrix(rnorm(n*p),nrow=n)
#'prob=1/exp(-(0.48+X%*% beta)+1) #ln(p/(1-p))=bX
#'y=1/(1+exp(-X%))
#'y=(prob>runif(n))
#'y= ifelse((prob>runif(n)), 1, 0)
#'data=cbind(y,X)
#' pqlLogist(data=data,G=G,nk=nk)
pqlLogist=function(data, G,nk){
  p=ncol(data[,-1]);n=nrow(data[,-1])
  y=data[,1];X=data[,-1]
  b <- runif(p, 0, 1)
  beta =matrix(b,nrow=p)
  Rm=matrix(rep(0, nk*G),ncol=G)
  mr=matrix(rep(0,G*nk), ncol=nk)
  R=matrix(rep(0,nk*n),ncol=n)
  be=matrix(rep(0,G*p),ncol=p);u=rep(0,G); v=rep(0,G);w=rep(0,G)
  for(i in 1:G ){
    mr[i,]=sample(1:n,nk,replace=F)
    r=matrix(c(1:nk,mr[i,],ncol=nk,byrow=T))
    R[t(r)]=1
    X0= R%*%X
    y0=R%*%y
    y0=as.factor(y0)
    be[i,]<-glm(y0~ X0-1, family = quasibinomial(link = logit))$coefficients
    u[i]=var(be[i,])
    v[i]=1/u[i]
    w[i]=(v[i])/(sum(v))}
  betaA=rep(0,p);betaW=rep(0,p)
  for (j in 1:p) {
    betaA[j]=sum(be[,j])/G
    betaW[j]=t(be[,j])%*%w}
  beta <- as.numeric(beta)
  betaA <- as.numeric(betaA)
  betaW <- as.numeric(betaW)
  var(beta-betaA);var(beta-betaW)
  return(list(betaW=betaW,betaA=betaA,MSEW=var(beta-betaW),MSEA=var(beta-betaA)))
}

