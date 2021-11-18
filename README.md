# QPLoF

This package carries out the **lack-of-fit test** for quantile regression processes. Quantile regression models are often used in practice, but model misspecifications may lead toincorrect statistical inferences. In our package, a lack-of-fit test for quantile regression processes is realized for those cases with multivariate covariates.

## Installation

You can use `devtools` to directly install the latest version from Github

```R
# install.packages("devtools")
devtools::install_github("lqiaochu/QPLoF")
```

## Example

Here is a toy example:

```r
## Parameter Settings
n = 200			# number of obsevations
p = 3			# number of covariates dimensions
l0 = 0.4			# start point of quantile interval
u0 = 0.6			# end point of quantile interval
tau_kn = 11			# number of grids for quantile interval
spline_kn = 5			# number of knots for B-splines
degree = 3			# order of B-splines
beta0 = c(1,2,1)			# initial coefficient vectors for beta
gamma0 = rep(1,p)			# initial coefficient vectors for gamma
B = 50			# bootstrap times (at least 1000 times, it's set to 50 just for testing)
v_index	= 1			# index of nonlinear functions

## Data generalization
epsilon = rnorm(n,0,1)			# drawn from standard normal distribution
# Here we use example functions in gen_data_example.R in our package to generate preditor X and response y.
# preditor X is drawn from exponential distribution exp(1)
# preditor y is generated by an exponential nonlinear function
X = gen_x(n = n, p = p)
y = gen_y(x = X, beta0 = beta0, gamma0 = gamma0, e = epsilon, v = v_index)			
tau = seq(l0,u0,length.out = tau_kn)			# quantile levels
breaks = seq(l0,u0,length.out = spline_kn)			# breaks for B-splines
```

Then we first estimate the quantile regression process based on the B-splines approximation,  we use the **MM** algorithm to obtain the solution of optimizaiton.  The result is stored in the variable `b_hat` ,  a vector correpond to the estimated coefficients. You can do the estimation by

```R
## prepare for approximation
Phi_tau = gen_bspline(tau = tau, breaks = breaks, degree = degree)
z_tau = kronecker(Phi_tau,X)
## optimization via MM algorithm
result = comp.B.MM.c(ztau = z_tau, Y = y, tau = tau, breaks = breaks, basis.order = degree, n = n, p = p, maxiter=200, tol=10^-8, epsilon=0.01) 
b_hat = as.vector(result$B)			# store the estimated coefficient matrix
e_hat = result$Residual			# store the estimated error

```

Next, you can easily calculate the **lack-of-fit** test statistic

```R
Phi = phi(y = Y, z_tau = z_tau, tau = tau, b = b_hat, n = n)
U_star = Un(x = X, phi = Phi)			# store the test statistic
```

where `Phi` is the parameter needed for calculating the test statistic.

To approximate the limiting distribution of the **lack-of-fit** statistic, we consider paired bootstrap to resample. You can calculate the  bootstrap statistic by

```R
## paired bootstrap
ztau,Y,tau,breaks,basis.order,n,p,
B_boot<-comp.B.MM.c(ztau = z_tau_boot, Y = Y_boot, tau = tau, breaks = breaks, basis.order = degree, n = n, p = p, maxiter=200, tol=10^-8, epsilon=0.01)$B
Phi_boot = phi(y = Y_boot, z_tau = z_tau_boot, tau = tau, b = as.vector(B_boot), n = n)
## calculate the boostrap test statistic
U_boot = Un.b.c(x = as.matrix(X), x.b = as.matrix(x_boot), phi = Phi, phi.b = Phi_boot)
```













