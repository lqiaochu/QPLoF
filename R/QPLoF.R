#' Title QPLoF
#'
#' @param X
#' @param y
#' @param tau
#' @param n
#' @param p
#' @param breaks
#' @param degree
#' @param maxiter
#' @param tol
#' @param epsilon
#' @param B
#'
#' @return pvalue
#' @export
#'
#' @examples
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

#
#
# test_sim<-function(seq)
# {
#   res = rep(0,length(seq))
#   for(j in 1:length(seq))
#   {
#     set.seed(seq[j])
#     ######## generate data ##########
#     epsilon = rnorm(n,0,1)
#     X = gen_x(n,p)
#     y = gen_y(X,beta0,gamma0,epsilon,v=v_index) #h0
#     Y = rep(y,length(tau))
#     Phi_tau = gen_bspline(tau,breaks,degree)
#     z_tau = kronecker(Phi_tau,X)
#     ###### estimation #######
#     result = comp.B.MM.c(z_tau,Y,tau,breaks,degree,n,p,maxiter=200,tol=10^-8,epsilon=0.01)
#     # result = comp.B.MM(z_tau,Y,n,p,tau,breaks,degree,maxiter=200,tol=10^-8,epsilon=0.01)
#     B_hat = result$B
#     b_hat = as.vector(B_hat)
#     e_hat = result$Residual
#     Phi = phi(Y,z_tau,tau,b_hat,n)
#     U_star = Un(X,Phi)
#     U_boot = vector(length = B)
#     for(i in 1:B)
#     {
#       print(i)
#       # U_boot[i] = Un_Boot(n,p,z_tau,X,b_hat,e_hat,breaks,degree,tau) ## wild bootstrap
#       # U_boot[i] = Un_Boot(n,p,y,X,breaks,degree,tau)
#       Ind<-sample(1:n,n,replace = T)
#       x_boot<-as.matrix(X)[Ind,]
#       y_boot<-y[Ind]
#       Y_boot = rep(y_boot,length(tau))
#       Phi_tau = gen_bspline(tau,breaks,degree)
#       z_tau_boot = kronecker(Phi_tau,x_boot)
#       B_boot<-comp.B.MM.c(z_tau_boot,Y_boot,tau,breaks,degree,n,p,maxiter=200,tol=10^-8,epsilon=0.01)$B
#       Phi_boot = phi(Y_boot,z_tau_boot,tau,as.vector(B_boot),n)
#       U_boot[i]<-Un.b.c(as.matrix(X),as.matrix(x_boot),Phi,Phi_boot)
#     }
#     res[j] = sum(U_star<U_boot)/B
#     # res is p_value
#   }
#   return(res)
# }


# ######## parameter setting ##########
# l0 = 0.4
# u0 = 0.6
# # spline_kn = 7 # knots for splines P=5,7
# spline_kn = 5 # knots for splines P=3
# tau_kn = 11 # knots for tau
# tau = seq(l0,u0,length.out = tau_kn)
# breaks<-seq(l0,u0,length.out = spline_kn)
# degree = 3
# n = 200
# p = 3
# v_index = 1
# # beta0 = c(1,1.5,1,2,1) #p=5
# # beta0 = c(2,1,4,2,1.5,2,1) #p=7
# beta0 = c(1,2,1) #p=3
# gamma0 = rep(1,p)
# B = 1000
# MC = 500
#
# ############plot for f1 and f2#############
# # p=3
# num = 50
# x1 = rep(1,num)
# x = seq(0,10,length.out = num)
# x3 = seq(0,10,length.out = num)
# f1 = x1^2+x^2+x3^2
# f2 = exp(0.1*(f1))
# par(pin = c(3.8,2.9))
# plot(1,1,xlim = c(0,10),ylim = c(0,500),type = 'n', ann = F, las=1)
# grid()
# lines(x,f1,lty=1,col = 'coral')
# lines(x,f2,lty=2,col = 'skyblue')
# # par(mai=c(2.5,0.5,0.5,2.5))
# legend(x = 'topleft', cex = 0.8,legend = c(expression(f[1]:f(x[i]) == x[i]^T*x[i]),expression(f[2]:f(x[i]) == exp(x[i]^T*x[i]/10))),lty = c(1,2),col = c('coral','skyblue'))
#


# #################估计画图###################
# est_beta = list()
# # real_betatau <-beta0[1]+gamma0[1]*qnorm(tau)
# # par(family='STKaiti')
# # plot(tau,real_betatau,type = "l",pch = 20,xaxt="n", yaxt="n",xlab = "分位数点",ylab = "参数估计值",lty = 1,col = 'black',lwd = 2)
# # axis(side=2,at=c(0.0,0.5,1.0,1.5,2.0))
# # axis(side=1,at=c(0.4,0.45,0.5,0.55,0.6))
# dyn.load("mybspline_c.so")
# for(t in 1:100)
# {
#   set.seed(t)
#   print(t)
#   epsilon = rnorm(n,0,1)
#   X = gen_x(n,p)
#   y = gen_y(X,beta0,gamma0,epsilon,v=0) #h0
#   Y = rep(y,length(tau))
#   Phi_tau = gen_bspline(tau,breaks,degree)
#   z_tau = kronecker(Phi_tau,X)
#   ###### estimation #######
#   result = comp.B.MM.c(z_tau,Y,tau,breaks,degree,n,p,maxiter=200,tol=10^-8,epsilon=0.01)
#   ############################# plot
#   colorvec = c("red","gold4","blue","green","orange","deeppink2","gray32" )
#   e.quant.p<-as.vector(result$B%*%t(Phi_tau))
#   e.quant.p = matrix(e.quant.p,nrow = ncol(X))
#   # lines(tau,e.quant.p[1,], type = "l", pch = 16,ylim = c(0,2),xlim = c(0.4,0.6),xlab = "分位数点",ylab = "参数估计值",lty = 3,col = colorvec[t%%7],lwd = 0.6)
#   est_beta[[t]] = e.quant.p
# }
# dyn.unload("mybspline_c.so")
# saveRDS(est_beta,file = "estbetan300p5.rds") # sin(2 x^2) setting
# ### n300p3 n300p5
#
# # plot
# estbeta = readRDS("/Users/mac/chuchu/quantile_process/code210729/result/estimateresult/estbetan300p5.rds")
# real_betatau <-beta0[1]+gamma0[1]*qnorm(tau)
# par(family='STKaiti')
# plot(tau,real_betatau,type = "l",pch = 20,ylim = c(1,3),xlab = "分位数点",ylab = "参数估计值",lty = 1,col = 'black',lwd = 2)
# axis(side=1,at=c(0.4,0.45,0.5,0.55,0.6))
# sum = vector(length = 11)
# for(i in 1:100)
# {
#   est_beta = estbeta[[i]][1,]
#   sum = sum+est_beta
#   colorvec = c("red","gold4","blue","green","orange","deeppink2","gray32" )
#   lines(tau,est_beta, type = "l", pch = 16,xlim = c(0.4,0.6),xlab = "分位数点",ylab = "参数估计值",lty = 3,col = colorvec[i%%7],lwd = 0.6)
# }
# mean_est_beta = sum/100
# mse = sum((mean_est_beta-real_betatau)^2)
# mse
#
# lines(tau,mean_est_beta,col="blue",lty = 2,type = "l",lwd = 2,pch = 17)
# # legend("topright", cex = 0.1, c("real","est"), pch=c(20, 16), col=c("black", "blue"),xjust = 0,yjust = 1,adj = 0,bty = "o")



# ########## 并行 ##########
# ###### our method ######
# no_cores <- 15
# cl <- makeCluster(no_cores)
# seq<-clusterSplit(cl,1:MC)
# # seq = 1:MC
# clusterExport(cl,ls())
# clusterCall(cl,function() library(fda))
# clusterCall(cl,function() library(combinat))
# clusterCall(cl,function() dyn.load("mybspline_c.so"))
#
# clusterExport(cl,varlist = c("l0","u0","tau","breaks","degree","n","p","beta0","gamma0","B"))
# starttime=Sys.time()
# par_res = parSapply(cl,seq,test_sim)
# saveRDS(par_res,file = paste("pairedbootexp01","n",n,"p",p,"h1.rds",sep = '-')) # exp(0.1 x^2) setting
# # saveRDS(par_res,file = paste("pairedbootx2","n",n,"p",p,"h0.rds",sep = '-'))
#
# endtime = Sys.time()
# clusterCall(cl,function() dyn.unload("mybspline_c.so"))
# stopCluster(cl)
#
# result = readRDS("/Users/mac/chuchu/quantile_process/code210729/newresult/pairedbootx2-n-200-p-3-h0.rds")
# size = sum(result<0.05)/500


# ############# ROC plot #################
# pvalues.nul = readRDS("/Users/mac/chuchu/quantile_process/code210729/result/pairedbootstrapresult/pairedbootx2-n-500-p-7-h0.rds")
# pvalues.alt.x2 = readRDS("/Users/mac/chuchu/quantile_process/code210729/result/pairedbootstrapresult/pairedbootexp01-n-500-p-7-h1.rds")
# pvalues.nul.esc = readRDS("/Users/mac/chuchu/quantile_process/code210729/result/JCE/JCEpairedbootx2-n-500-p-7-h0.rds")
# pvalues.alt.x2.esc = readRDS("/Users/mac/chuchu/quantile_process/code210729/result/JCE/JCEpairedbootexp01-n-500-p-7-h1.rds")
#
# p.value=function(alpha,value)
# {
#   sum(value<alpha)/length(value)
# }
# sep=as.matrix(seq(0,1,by=0.01))
# x0=apply(sep,1,p.value,value=pvalues.nul)
# y0=apply(sep,1,p.value,value=pvalues.alt.x2)
#
# x1=apply(sep,1,p.value,value=pvalues.nul.esc)
# y1=apply(sep,1,p.value,value=pvalues.alt.x2.esc)
#
# # tit = substitute(tau==a, list(a = tau))
# par(pin = c(4,2))
# plot(cbind(x0,y0),type='n', xaxt = "n",yaxt = "n", xlab='False positive rate',ylab='True positive rate')
# grid()
# axis(side=1, at = seq(0,1,0.25))
# axis(side=2, at = seq(0,1,0.25))
# lines(cbind(x0,y0),col=1,lty=1,lwd=2)
# lines(cbind(x1,y1),col=1,lty=2,lwd=2)
# #abline(0,1,lty=1,col=3,lwd=1)
# # legend("bottomright",legend=c('Proposed','QRS','Boot'),lty=c(1,2,3),lwd=2,seg.len=5)














