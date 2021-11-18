#' Title comp.B.MM.c
#' compute the B-spline coef matrix B using MM algorithm via C code
#' @param ztau transformed preditor variable
#' @param Y response variable
#' @param tau quantile levels
#' @param breaks breaks for B-splines, including starting and finishing points of the interval
#' @param basis.order order of B-splines
#' @param n number of observations
#' @param p dimension of covariate
#' @param maxiter maximum number of iteration during MM algorithm
#' @param tol tolerance for convergence
#' @param epsilon parameter used in MM, the default value is 0.001
#'
#' @return
#' @export list
#'
#' @examples
comp.B.MM.c<-function(ztau,Y,tau,breaks,basis.order,n,p,maxiter=1000,tol=10^-8,epsilon=0.001)
{
  kn<-length(tau)-1
  m<-length(breaks)+basis.order-2   # m: number of basis functions

  ### prepare the data to get the standard form for MM algorithm ###
  theta<-solve(t(ztau)%*%ztau)%*%t(ztau)%*%Y   # theta: strating point of theta
  ### the MM algorithm ###
  theta<-.Call("compbmmc",ztau,Y,rep(tau,each=n),theta,as.integer(n),as.integer(p*m),as.integer(kn+1),as.integer(maxiter),tol,epsilon)
  R<-Y-ztau%*%theta

  B<-matrix(theta,nrow=p)
  R<-matrix(R,ncol=kn+1)
  res<-list(B=B,Residual=R)
  return(res)
}
