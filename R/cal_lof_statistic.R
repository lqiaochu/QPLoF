
#' Title Un
#' calculate the lack-of-fit test statistic
#' @param x predictor variable
#' @param phi results of funciton Phi with respect to x, in vector form
#'
#' @return number
#' @export
#'
#' @examples
Un<-function(x,phi)
{
  n = nrow(x)
  ind = combn(n,2)
  res<-apply(ind,2,function(ind,x,phi) {
    u<-prod(phi[ind])
    v<-(x[ind[1],]-x[ind[2],])^2
    v<-sqrt(sum(v))
    return(u*v)
  },x,phi)
  res = -sum(res)/ncol(ind)
  return(res)
}


#' Title Un.b.c
#' calculate the bootstrap lack-of-fit test statistic
#' @param x original sample
#' @param x.b boostrap sample of x
#' @param phi results of funciton Phi with respect to x, in vector form
#' @param phi.b results of funciton Phi with respect to x.b, in vector form
#'
#' @return number
#' @export
#'
#' @examples
Un.b.c<-function(x,x.b,phi,phi.b)   # Bootstrap LOF test statistic
{
  # in this function, we use ||x_j-x_i||=sqrt(sum((x_j-x_i)^2))
  # in this function, we use c code to compute the result
  x<-as.matrix(x)
  x.b<-as.matrix(x.b)
  n<-nrow(x)
  p<-ncol(x)
  res<-.Call("unbc",x,x.b,phi,phi.b,as.integer(n),as.integer(p))
  return(res)
}
