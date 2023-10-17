source("simulation/generate_data.R")
n <- 50
p <- 100
beta_star <- c(1.2, -0.8, 0.7, 0, 0, -1.5, -1, 1.4, rep(0, p - 8))
theta <- c(rep(log(0.5 * p), 5), rep(0, p - 5))
sigma <- 0.5
rho <- 0.5


# data <- generate_multinom_data(n, p, beta_star, sigma, rho, theta = theta)
# sig <- MC_varB(n, p, beta_star, sigma, rho, theta = theta, N_MK = 1000000 %/% p, type = "multinom")
# sqrt(sig[c(1, 6), c(1, 6)])

# W <- AR(n, p, rho) + rep(theta, each = n)
# X <- exp(W) / rowSums(exp(W))
# X[1, c(1, 50)]


MC_gen_varB <- function(n, p, beta_star, sigma, rho, theta = NA, N_MK = 10000 %/% p) {
    varB <- array(NA, c(N_MK, p, p))
    for (mk in 1:N_MK) {
        tmp_data <- generate_dirmult_data(n, p, beta_star, sigma, rho, theta = theta)
        B <- tmp_data$Z - tmp_data$X
        varB[mk,,] <- var(B)
    }
    Sig_B <- apply(varB, c(2, 3), mean)

    return(Sig_B)
}

generate_dirmult_data <- function(n, p, beta_star, sigma, rho, theta = theta) {
    W <- AR(n, p, rho) + rep(theta, each = n)
    X <- exp(W) / rowSums(exp(W))
    for (i in 1:n) {
        Q <- rdirichlet(1, alpha = 5e+03 * X[i, ])
        nSeq <- rnbinom(n, prob = 0.01, size = 300 / 0.99)
        W[i, ] <- rmultinom(1, nSeq, prob = Q)[, 1]
    }
    W <- W + 0.5 ## recommendated by 2023 BKA paper
    Z <- W / rowSums(W)
    Z <- log(Z)
    X <- log(X)
    y <- X %*% beta_star + rnorm(n, sd = sigma)
    return(list(X = X, Z = Z, y = y))
}

data <- generate_dirmult_data(n, p, beta_star, sigma, rho, theta = theta)
sig <- MC_gen_varB(n, p, beta_star, sigma, rho, theta = theta, N_MK = 10000 %/% p)
sqrt(abs(sig[1:6,1:6]))