#############################################################
####### code for high-dimensional calibration model
#############################################################
source('simulation/functions.R')
source('simulation/generate_data.R')
# 1. generate data--------
n = 200 # sample size
p = 300 # dimension
## parameters values
s2_u = 0.5
mu_u = rep(-0.5*s2_u,p) # mean of log u
Sigma_u = s2_u*diag(rep(1,p)) # variance of log u
mu_x = round(c(rep(log(p/2),5),rep(0,p-5)),2) # mean of log x
Sigma_x = matrix(NA,p,p) # variance of log x
rho = 0.2
for(j in 1:p){
  for(l in 1:p){
    Sigma_x[j,l] = rho^(abs(j-l))
  }
}
mu_w = mu_x + mu_u # mean of log w
Sigma_w = Sigma_x + Sigma_u # variance of log w
alpha_real = rep(0,p-1) # true coefficient, (p-1)-dim 
alpha_real[1:7] = c(1, -0.8, 1.5, 0.6, -0.9, 1.2, 0.4) 
## generate dataset
set.seed(123)
X = exp(mvrnorm(n,mu_x,Sigma_x))  # n * p
U = exp(mvrnorm(n,mu_u,Sigma_u)) 
W = X*U 
Z = X # true compositional data
V = W # observed compositional data
for(i in 1:n){
  Z[i,] = Z[i,]/sum(X[i,])
  V[i,] = V[i,]/sum(W[i,])
}
## covariates in log-contrast model
ZZ = VV = matrix(NA,n,p-1) 
for(j in 1:(p-1)){
  ZZ[,j] = log(Z[,j]/Z[,p]) # true covariates
  VV[,j] = log(V[,j]/V[,p]) # observed covariates
}
## generate response
eps1 = rnorm(n,0,0.5) 
y = as.vector(ZZ%*%alpha_real) + eps1 # response,  n-dim

# 2. estimate alpha_j using our proposed method------------
## input
VV = VV # observed covariates in log-contrast model with measurement error, a n*(p-1) matrix
y = y # response, a n-dim vector
EstimateSigma = F # TRUE: estimated Sigma_x and mu_x; FALSE: known Sigma_x and mu_x
Shrinkage = T # TRUE: shrinkage estimator, FALSE: node-wise estimator
## run our proposed method
fit_prop = fn_proposed_original(VV,y,EstimateSigma,Shrinkage)$beta
fit_lasso = fn_lasso(VV,y,EstimateSigma,Shrinkage)$beta
fit_delasso = fn_delasso(VV,y,EstimateSigma,Shrinkage)$beta

trunc = function(l){
  l[1:10]
}
print(trunc(alpha_real))
print(trunc(fit_prop))
print(trunc(fit_lasso))
print(trunc(fit_delasso))
print(eval_results(fit_prop, ZZ, alpha_real))
print(eval_results(fit_lasso, ZZ, alpha_real))
print(eval_results(fit_delasso, ZZ, alpha_real))
## output
## the true parameter, point estimator, variance estimator and p-value for alpha_j, j = 1,...,p-1.









