source("simulation/generate_data.R")
n = 100
p = 800
# simulation settings
data_type <- "lognormal"
# data_type <- "dirichlet"
# data_type <- "multinom"
# data_type <- "dirmult"
# create a list of different n and p values


sigma <- 0.5
rho <- 0.5
tau <- 0.5
# tau_list <- c(10.5,11.5,12.5)

beta_star <- c(1.2, -0.8, 0.7, 0, 0, -1.5, -1, 1.4, rep(0, p - 8))
# theta <- c(rep(log(0.01 * p), 100), rep(0, p - 100))
theta <- c(rep(log(0.2 * p), 5), rep(0, p - 5))

data <- generate_data(n, p, beta_star, sigma, tau, rho, theta, type = data_type)

sigx = cov(data$X)
dim(sigx)
source("R/ADMM_proj.R")
library(Rcpp)
sourceCpp("R/ADMM_proj.cpp")

proj.c = ADMM_proj_cpp(mat=sigx)$mat
proj.r = ADMM_proj(mat=sigx)$mat
norm(proj.c-proj.r, 'F')




#################### lasso covariance con  ##########
library(Rcpp)
sourceCpp("/storage/home/wkt5100/work/Eric-Lasso/R/lasso_covariance_constrained.cpp")

# Example usage

n <- 50
p <- 100
X <- matrix(rnorm(n * p), n, p)
beta_star <- c(1.2, -0.8, 0.7, 0, 0, -1.5, -1, 1.4, rep(0, p - 8))
y <- X %*% beta_star + rnorm(n)
XX <- t(X) %*% X / n
Xy <- t(X) %*% y / n
lambda <- 0.5


beta_start <- rep(0, p)

library(Rcpp)
sourceCpp("R/lasso_covariance.cpp")

coef.c = lasso_covariance_rcpp(n,p, 0.01, XX, Xy, beta_start)$coefficients
eval_results(coef.c, X, beta_star)

source("R/lasso_covariance.R")

coef.r = lasso_covariance(n,p,lambda= 0.01, XX=XX,Xy= Xy, beta.start = beta_start, penalty='lasso')$coefficients
eval_results(coef.r, X, beta_star)

# cl <- makeCluster(5)
# clusterExport(cl, c("n", "p", "lambda", "control", "XX", "Xy", "beta_start", "X","y"))
# clusterEvalQ(cl, library(Rcpp))
# clusterEvalQ(cl, sourceCpp("/storage/home/wkt5100/work/Eric-Lasso/R/lasso_covariance_constrained.cpp"))
# clusterExport(cl, c("lasso_covariance_con_cpp"))

# source("R/lasso_covariance_constrained.R")
# clusterExport(cl, c("lasso_covariance_con","soft_threshold"))
# result.c <- parLapply(cl, 1:5, function(i) lasso_covariance_con(n, p, lambda, control, XX, Xy, beta_start))
# # result.c <- parLapply(cl, 1:5, function(i) lasso_covariance_con_cpp(n, p, lambda, control, XX, Xy, beta_start))
# stopCluster(cl)
# coef.c = result.c[[1]]$coefficients
# eval_results(coef.c, X, beta_star)
# coef.c = result.c[[2]]$coefficients
# eval_results(coef.c, X, beta_star)






# result.c <- lasso_covariance_con_cpp(n, p, lambda, control, XX, Xy, beta_start)
# coef.c = result.c$coefficients
# eval_results(coef.c, X, beta_star)

# source("R/lasso_covariance_constrained.R")

# result.r <- lasso_covariance_con(n, p, lambda, control, XX, Xy, beta_start)
# coef.r = result.r$coefficients
# eval_results(coef.r, X, beta_star)