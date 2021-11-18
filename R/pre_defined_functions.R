#' Title phi
#'
#' @param y response variable
#' @param z_tau transformed preditor variable
#' @param tau quantile levels
#' @param b coefficient value
#' @param n number of observations
#'
#' @return number
#' @export
#'
#' @examples
phi<-function(y,z_tau,tau,b,n)
{
  kn = length(tau)-1
  Tau = matrix(rep(tau,each=n),nrow = n)
  psi = y - z_tau%*%b
  psi = matrix(psi,nrow = n)
  psi = Tau - (psi<0)
  return(rowSums(psi))
}


#' Title gen_bspline
#' generate b-splines based on tau
#' @param tau  quantile level seq
#' @param breaks  breaks for B-splines
#' @param degree  basis order of B-splines
#'
#' @return vector
#' @export
#'
#' @examples
gen_bspline<-function(tau,breaks,degree)
{
  # splines_tau = bs(tau,knots = breaks,degree = degree)
  splines_tau = bsplineS(tau,breaks,degree)
  return(splines_tau)
}
