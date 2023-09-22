library(dirmult)

n=50
p=50
rho=0.5
tau=0.5



AR = function(n,p,rho)
{
  z=matrix(NA,n,p)
  for (j in 1:(p))
  {
    z[,j]=rnorm(n)
 }
X=matrix(NA,n,p)
  X[,1]=z[,1]
  for (j in 2:(p)){
    X[,j]=rho*X[,j-1]+sqrt(1-rho^2)*z[,j]
  }  
  return(X)
}

get_std <- function(N_sim){
  varB <- matrix(NA,N_sim,p)
  for(N in 1:N_sim){
    W1 <-AR(n,p,rho)
    u <- c(rep(log(0.5*p),5),rep(0,p-5))
    W <-  W1 + rep(u,each=n)
    X <- matrix(NA,n,p)
    for (i in 1:n){
      for (j in 1:p){
        X[i,j] = exp(W[i,j])/sum(exp(W[i,]))
      }}
    # rowSums(X)
    W=X ## Count matrix simulated from community compositions
    for(i in 1:nrow(W)){
      nSeq <- rnbinom(1, mu = 10000, size = 25)
      W[i, ] <- rmultinom(1, nSeq, prob=X[i, ])[, 1]
    }
    W=W+0.5   ## recommendated by 2023 BKA paper
    Z=W
    for(i in 1:nrow(W)){
      Z[i,]=W[i,]/rowSums(W)[i]
    }
    # rowSums(Z)
    Z <- log(Z)
    X <- log(X)
    B <- Z-X
    for(i in 1:ncol(B)){
      varB[N,i]=sd(B[,i])
    }
  }
  colMeans(varB)
}

std <- get_std(10000)
mean(std[1:5])
mean(std[6:p])


## Case 1: Using Logistic-Normal distribution to simulate heterogeneous communities
# W1 <-AR(n,p,rho)
# u <- c(rep(log(0.5*p),5),rep(0,p-5))
# W <-  W1 + rep(u,each=n)
# X <- matrix(NA,n,p)
# for (i in 1:n){
# for (j in 1:p){
# 	X[i,j] = exp(W[i,j])/sum(exp(W[i,]))
# }}
# rowSums(X)
# B <- matrix(rnorm(n * p, mean = 0, sd = tau), nrow = n, ncol = p)
# Sig_B <- tau**2 * diag(p)
# Z <- X * exp(B)
# Z <- Z / rowSums(Z)
# rowSums(Z)
# Z <- log(Z)
# X <- log(X)


# ## Case 2: Using Dirichlet distribution to simulate homogeneous communities
# X=matrix(NA,n,p)
# for(i in 1:n){X[i,]=rdirichlet(1,alpha=rep(1/p,p))}
# rowSums(X)
# B <- matrix(rnorm(n * p, mean = 0, sd = tau), nrow = n, ncol = p)
# Sig_B <- tau**2 * diag(p)
# Z <- X * exp(B)
# Z <- Z / rowSums(Z)
# rowSums(Z)
# Z <- log(Z)
# X <- log(X)


# ## Case 3: Simulation data agaist model assumptions Z=X+B
# ## Using Multinomial to simulate observed counts
# W1 <-AR(n,p,rho)
# u <- c(rep(log(0.5*p),5),rep(0,p-5))
# W <-  W1 + rep(u,each=n)
# X <- matrix(NA,n,p)
# for (i in 1:n){
# for (j in 1:p){
# 	X[i,j] = exp(W[i,j])/sum(exp(W[i,]))
# }}
# rowSums(X)
# W=X ## Count matrix simulated from community compositions
# for(i in 1:nrow(W)){
# nSeq <- rnbinom(1, mu = 10000, size = 25)
# W[i, ] <- rmultinom(1, nSeq, prob=X[i, ])[, 1]
# }
# W=W+0.5   ## recommendated by 2023 BKA paper
# Z=W
# for(i in 1:nrow(W)){
# Z[i,]=W[i,]/rowSums(W)[i]
# }
# rowSums(Z)
# Z <- log(Z)
# X <- log(X)



# ## Case 4: Simulation data agaist model assumptions Z=X+B
# ## Using Multinomial to simulate observed counts
# X=matrix(NA,n,p)
# for(i in 1:n){X[i,]=rdirichlet(1,alpha=rep(1/p,p))}
# rowSums(X)
# W=X ## Count matrix simulated from community compositions
# for(i in 1:nrow(W)){
# nSeq <- rnbinom(1, mu = 10000, size = 25)
# W[i, ] <- rmultinom(1, nSeq, prob=X[i, ])[, 1]
# }
# W=W+0.5   ## recommendated by 2023 BKA paper
# Z=W
# for(i in 1:nrow(W)){
# Z[i,]=W[i,]/rowSums(W)[i]
# }
# rowSums(Z)
# Z <- log(Z)
# X <- log(X)



