#############################################################
####### functions for high-dimensional calibration model
#############################################################
library(MASS)
library(glmnet)
library(corpcor)
library(lars)
library(scalreg)
## Node-wise estimator of covariance matrix based on javanmard2018
fn_cov <- function(X){
  for(j in 1:p){  X[,j] = X[,j] - mean(X[,j]) }
  C = matrix(1,p,p)
  lambda = tau2 = rep(NA,p)
  for(j in 1:p){
    lambdas <- seq(0,1, length.out = 1000)
    lasso_model <- cv.glmnet(X[,-j], X[,j], alpha = 1, nfolds = 10, family = "gaussian")
    lambda[j] = lasso_model$lambda.min
    lasso_best <- glmnet(X[,-j], X[,j], alpha = 1, lambda = lambda[j], family = "gaussian")
    coef = as.vector(lasso_best$beta)
    if(length(coef)==p){# lasso omit intercept if predictor and response are zero-mean
      C[j,-j] = -coef[-1]
      resid = X[,j] - coef[1] - X[,-j]%*%coef[-1]
    }else{
      C[j,-j] = -coef
      resid = X[,j] - X[,-j]%*%coef
    }
    tau2[j] = t(resid)%*%X[,j]/n
  }
  T2 = diag(tau2)
  ## estimated covariance matrix
  out = solve(C)%*%T2
  return(out)
}
## function for the proposed method
fn_proposed <- function(VV,y, alpha_real, W, mu_u, Sigma_u, EstimateSigma=T,Shrinkage=T, Noestimate.data=NA, p_val = 0.05){
  n = dim(VV)[1]
  p = dim(VV)[2] + 1
  mean.Z = rep(0, p)
  ## obtain the mu_x and Sigma_x
  if(EstimateSigma){
    ## estimate mu_x, mu_w, Sigma_x, Sigma_w 
    ## assume s2_u is known, that is, mu_u and Sigma_u are known
    mu_w_hat = apply(log(W),2,mean)
    mu_x_hat = mu_w_hat - mu_u
    if(Shrinkage){
      Sigma_w_hat = cov.shrink(log(W)) # p * p
    }else{
      Sigma_w_hat = fn_cov(log(W)) # p * p
    }
    Sigma_x_hat = Sigma_w_hat - Sigma_u
    ## use estimated nuisance parameters
    mu_x = mu_x_hat; mu_w = mu_w_hat
    Sigma_x = Sigma_x_hat; Sigma_w = Sigma_w_hat
  }else {
    mu_x = Noestimate.data$mu_x; mu_w = Noestimate.data$mu_w
    Sigma_x = Noestimate.data$Sigma_x; Sigma_w = Noestimate.data$Sigma_w
  }
  ## obtain mu(V) and Sigma
  mu_1 = mu_x[1:(p-1)]-mu_x[p] # mean of tilde(X-p) or Z-p (p-1) * 1 
  mu_2 = mu_w[1:(p-1)]-mu_w[p] # mean of tilde(W-p) or V-p (p-1) * 1
  Sigma_11 = matrix(NA,p-1,p-1) # variance of tilde(X-p) or Z-p
  Sigma_22 = matrix(NA,p-1,p-1) # variance of tilde(W-p) or V-p
  for(j in 1:(p-1)){
    for(l in 1:(p-1)){
      Sigma_11[j,l] = Sigma_x[j,l] - Sigma_x[j,p] - Sigma_x[p,l] + Sigma_x[p,p]
      Sigma_22[j,l] = Sigma_w[j,l] - Sigma_w[j,p] - Sigma_w[p,l] + Sigma_w[p,p]
    }
  }
  B = Sigma_11%*%solve(Sigma_22) # "projection" Sigma_12=Sigma_11 (p-1) * (p-1)
  ZV = matrix(NA,n,p-1) # mu(V)
  for(i in 1:n){
    ZV[i,] = mu_1 + as.vector(B%*%(log(W[i,-p])-log(W[i,p])-mu_2))
  }
  Sigma_zv = Sigma_11%*%solve(Sigma_22)%*%Sigma_11 # Var(mu(V)) (p-1) * (p-1)
  ## centered dataset
  y = y - mean(y)
  for(j in 1:(p-1)){
    ZV[,j] = ZV[,j] - mean(ZV[,j]) 
  }
  ## calibration model
  rclasso_model <- cv.glmnet(ZV, y, alpha = 1, nfolds = 5, family = "gaussian", intercept = F)
  rclasso_best <- glmnet(ZV, y, alpha = 1, lambda = rclasso_model$lambda.min, family = "gaussian",intercept = F)
  alpha_rclasso = as.vector(rclasso_best$beta)
  ## proposed method
  # M = solve(Sigma_zv)
  M = solve(Sigma_zv + 1e-4 * diag(p-1))
  alpha_drclasso = alpha_rclasso + 1/n*M%*%t(ZV)%*%(y - ZV%*%alpha_rclasso)
  fit = scalreg(ZV,y) # scaled lasso estimator for sigma
  s2_eps_drclasso_hat = (fit$hsigma)^2
  Sigma_drclasso = (s2_eps_drclasso_hat/n)*M%*%(t(ZV)%*%ZV/n)%*%M
  se_drclasso = sqrt(diag(Sigma_drclasso)) 
  p_drclasso = 2*pnorm(abs(alpha_drclasso)/se_drclasso, lower.tail = F) # 2*P(T>t|H_0)
  ## output
  ### estimate, se, p.value (H_0: alpha_j = 0)
  # result = round(cbind(alpha_real, alpha_drclasso, se_drclasso, p_drclasso),3)
  # colnames(result) = c("real", "estimate", "se", "p.value")
  # rownames(result) = paste0("alpha_",1:(p-1))
  beta.opt=c(alpha_drclasso, -sum(alpha_drclasso))
  alpha_test = alpha_drclasso
  alpha_test[p_drclasso >= p_val] = 0
  tmp = -sum(alpha_test)
  alpha_test = c(alpha_test, tmp)
  result <- list(
    lambda.opt = rclasso_model$lambda.min,
    beta.opt = beta.opt,
    beta.lasso = c(alpha_rclasso, -sum(alpha_rclasso)),
    beta.test = alpha_test,
    p.beta = p_drclasso,
    mean.Z = mu_w,
    mean.y = mean(y)
  )
  return(result)
}

fn_proposed_original <- function(VV,y,EstimateSigma=F,Shrinkage=T){
  n = dim(VV)[1]
  p = dim(VV)[2] + 1
  ## obtain the mu_x and Sigma_x
  if(EstimateSigma){
    ## estimate mu_x, mu_w, Sigma_x, Sigma_w 
    ## assume s2_u is known, that is, mu_u and Sigma_u are known
    mu_w_hat = apply(log(W),2,mean)
    mu_x_hat = mu_w_hat - mu_u
    if(Shrinkage){
      Sigma_w_hat = cov.shrink(log(W)) # p * p
    }else{
      Sigma_w_hat = fn_cov(log(W)) # p * p
    }
    Sigma_x_hat = Sigma_w_hat - Sigma_u
    ## use estimated nuisance parameters
    mu_x = mu_x_hat; mu_w = mu_w_hat
    Sigma_x = Sigma_x_hat; Sigma_w = Sigma_w_hat
  }
  ## obtain mu(V) and Sigma
  mu_1 = mu_x[1:(p-1)]-mu_x[p] # mean of tilde(X-p) or Z-p (p-1) * 1 
  mu_2 = mu_w[1:(p-1)]-mu_w[p] # mean of tilde(W-p) or V-p (p-1) * 1
  Sigma_11 = matrix(NA,p-1,p-1) # variance of tilde(X-p) or Z-p
  Sigma_22 = matrix(NA,p-1,p-1) # variance of tilde(W-p) or V-p
  for(j in 1:(p-1)){
    for(l in 1:(p-1)){
      Sigma_11[j,l] = Sigma_x[j,l] - Sigma_x[j,p] - Sigma_x[p,l] + Sigma_x[p,p]
      Sigma_22[j,l] = Sigma_w[j,l] - Sigma_w[j,p] - Sigma_w[p,l] + Sigma_w[p,p]
    }
  }
  B = Sigma_11%*%solve(Sigma_22) # "projection" Sigma_12=Sigma_11 (p-1) * (p-1)
  ZV = matrix(NA,n,p-1) # mu(V)
  for(i in 1:n){
    ZV[i,] = mu_1 + as.vector(B%*%(log(W[i,-p])-log(W[i,p])-mu_2))
  }
  Sigma_zv = Sigma_11%*%solve(Sigma_22)%*%Sigma_11 # Var(mu(V)) (p-1) * (p-1)
  ## centered dataset
  y = y - mean(y)
  for(j in 1:(p-1)){
    ZV[,j] = ZV[,j] - mean(ZV[,j]) 
  }
  ## calibration model
  rclasso_model <- cv.glmnet(ZV, y, alpha = 1, nfolds = 5, family = "gaussian", intercept = F)
  rclasso_best <- glmnet(ZV, y, alpha = 1, lambda = rclasso_model$lambda.min, family = "gaussian",intercept = F)
  alpha_rclasso = as.vector(rclasso_best$beta)
  ## proposed method
  M = solve(Sigma_zv)
  alpha_drclasso = alpha_rclasso + 1/n*M%*%t(ZV)%*%(y - ZV%*%alpha_rclasso)
  fit = scalreg(ZV,y) # scaled lasso estimator for sigma
  s2_eps_drclasso_hat = (fit$hsigma)^2
  Sigma_drclasso = (s2_eps_drclasso_hat/n)*M%*%(t(ZV)%*%ZV/n)%*%M
  se_drclasso = sqrt(diag(Sigma_drclasso)) 
  p_drclasso = 2*pnorm(abs(alpha_drclasso)/se_drclasso, lower.tail = F) # 2*P(T>t|H_0)
  ## output
  ### estimate, se, p.value (H_0: alpha_j = 0)
  result = round(cbind(alpha_real, alpha_drclasso, se_drclasso, p_drclasso),3)
  colnames(result) = c("real", "estimate", "se", "p.value")
  rownames(result) = paste0("alpha_",1:(p-1))
  fit <- list(
    lambda.opt = rclasso_model$lambda.min,
    beta = alpha_drclasso
  )
  # return(result)
  return(fit)
}


fn_lasso <- function(VV,y,EstimateSigma=F,Shrinkage=T){
  n = dim(VV)[1]
  p = dim(VV)[2] + 1
  ## obtain the mu_x and Sigma_x
  rclasso_model <- cv.glmnet(VV, y, alpha = 1, nfolds = 5, family = "gaussian", intercept = F)
  rclasso_best <- glmnet(VV, y, alpha = 1, lambda = rclasso_model$lambda.min, family = "gaussian",intercept = F)
  alpha_rclasso = as.vector(rclasso_best$beta)
  ## proposed method
  ## output
  ### estimate, se, p.value (H_0: alpha_j = 0)
  result = round(cbind(alpha_real, alpha_rclasso),3)
  colnames(result) = c("real", "estimate")
  rownames(result) = paste0("alpha_",1:(p-1))
    fit <- list(
    lambda.opt = rclasso_model$lambda.min,
    beta = alpha_rclasso
  )
  # return(result)
  return(fit)
}
fn_delasso <- function(VV,y,EstimateSigma=F,Shrinkage=T){
  n = dim(VV)[1]
  p = dim(VV)[2] + 1
  ## obtain the mu_x and Sigma_x
  rclasso_model <- cv.glmnet(VV, y, alpha = 1, nfolds = 5, family = "gaussian", intercept = F)
  rclasso_best <- glmnet(VV, y, alpha = 1, lambda = rclasso_model$lambda.min, family = "gaussian",intercept = F)
  alpha_rclasso = as.vector(rclasso_best$beta)
  ## proposed method
  ## output
  ### estimate, se, p.value (H_0: alpha_j = 0)
  Sigma_22 = matrix(NA,p-1,p-1) # variance of tilde(W-p) or V-p
  for(j in 1:(p-1)){
    for(l in 1:(p-1)){
      Sigma_22[j,l] = Sigma_w[j,l] - Sigma_w[j,p] - Sigma_w[p,l] + Sigma_w[p,p]
    }
  }
  M = solve(Sigma_22)
  alpha_drclasso = alpha_rclasso + 1/n*M%*%t(VV)%*%(y - VV%*%alpha_rclasso)
  result = round(cbind(alpha_real, alpha_rclasso),3)
  colnames(result) = c("real", "estimate")
  rownames(result) = paste0("alpha_",1:(p-1))
    fit <- list(
    lambda.opt = rclasso_model$lambda.min,
    beta = alpha_drclasso
  )
  # return(result)
  return(fit)
}