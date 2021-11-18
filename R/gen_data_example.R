#' Title gen_x
#' example for generate preditor variable (drawn from exponential distribution exp(1))
#' @param n number of observations
#' @param p dimensions of covariate
#'
#' @return matrix
#' @export
#'
#' @examples
gen_x<-function(n,p)
{
  x1 = rep(1,n)
  for (i in 1:(p-1)){
    x = rexp(n,rate = 1)
    x1 = rbind(x1,x)
  }
  return(t(x1))
}


#' Title gen_y
#' example for generate predictor y (a nonlinear or linear form)
#' @param x preditor variable
#' @param beta0 initial value of coefficient beta
#' @param gamma0 initial value of coefficient gamma
#' @param e epsilon(error term)
#' @param v index for nonlinear part
#'
#' @return vector
#' @export
#'
#' @examples
gen_y<-function(x,beta0,gamma0,e,v)
{
  y = x%*%beta0+v*exp(0.1*rowSums(x^2))+(x%*%gamma0)*e  #setting2 expx2
  return(y)
}
