#' Title QPLoF
#' calculate pvalue of our lack-of-fit test based on paired bootstrap
#' @param X preditor variable
#' @param y response variable
#' @param tau quantile levels
#' @param n number of observations
#' @param p dimensions of covariate
#' @param breaks breaks for B-splines, including starting and finishing points of the interval
#' @param degree order of B-splines
#' @param maxiter maximum number of iteration during MM algorithm
#' @param tol tolerance for convergence
#' @param epsilon parameter used in MM, the default value is 0.001
#' @param B bootstrap times
#'
#' @return pvalue
#' @export
#'
#' @examples pavlue = QPLoF(X = X, y = y, tau = tau ,n = n , p = p, breaks = breaks, degree = degree, maxiter=1000, tol=10^-8, epsilon=0.001)
QPLoF<-function(X, y, tau ,n , p, breaks, degree, maxiter, tol, epsilon, B)
{
  Y = rep(y,length(tau))
  Phi_tau = gen_bspline(tau,breaks,degree)
  z_tau = kronecker(Phi_tau,X)

  # estimation
  result = comp.B.MM.c(z_tau,Y,tau,breaks,degree,n,p,maxiter=200,tol=10^-8,epsilon=0.01)
  B_hat = result$B
  b_hat = as.vector(B_hat)
  e_hat = result$Residual

  # calculate lack-of-fit test
  Phi = phi(Y,z_tau,tau,b_hat,n)
  U_star = Un(X,Phi)

  # paired bootstrap to resample
  U_boot = vector(length = B)
  for(i in 1:B)
  {
    Ind<-sample(1:n,n,replace = T)
    x_boot<-as.matrix(X)[Ind,]
    y_boot<-y[Ind]
    Y_boot = rep(y_boot,length(tau))
    Phi_tau = gen_bspline(tau,breaks,degree)
    z_tau_boot = kronecker(Phi_tau,x_boot)
    B_boot<-comp.B.MM.c(z_tau_boot,Y_boot,tau,breaks,degree,n,p,maxiter=200,tol=10^-8,epsilon=0.01)$B
    Phi_boot = phi(Y_boot,z_tau_boot,tau,as.vector(B_boot),n)
    U_boot[i]<-Un.b.c(as.matrix(X),as.matrix(x_boot),Phi,Phi_boot)
  }
  pvalue = sum(U_star<U_boot)/B
  return(pvalue)
}














