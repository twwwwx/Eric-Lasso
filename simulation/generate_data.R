library(MASS)
library(dirmult)
#--------------------------------------------
# Generate data
#--------------------------------------------

generate_data <- function(n, p, beta_star, sigma, tau, rho, theta = NA, cov_type = "AR", type = "lognormal") {
    # n: number of rows
    # p: number of columns
    # beta_star: true beta
    # sigma: noise level
    # tau: sd of measurement eror
    # rho: correlation between columns of W
    # theta: mean of W

    if (anyNA(theta)) theta <- c(rep(log(0.5 * p), 5), rep(0, p - 5))

    # Generate data
    if (cov_type == "AR") {
        W <- AR(n, p, rho) + rep(theta, each = n)
    } else if (cov_type == "CS") {
        Sigma <- sym_covariance_matrix(p, rho)
        W <- mvrnorm(n, mu = theta, Sigma = Sigma)
    } else {
        stop("cov_type must be either AR or CS")
    }
    B <- matrix(rnorm(n * p, mean = 0, sd = tau), nrow = n, ncol = p)
    Sig_B <- tau**2 * diag(p)
    if (type == "lognormal") {
        X <- exp(W) / rowSums(exp(W))
    } else if (type == "dirichlet") {
        X <- matrix(NA, n, p)
        for (i in 1:n) {
            X[i, ] <- rdirichlet(1, alpha = rep(1 / p, p))
            X[i, ][X[i, ] < 1e-10] <- 1e-10
        }
    } else {
        stop("type must be either lognormal or dirichlet")
    }
    Z <- X * exp(B)
    Z <- Z / rowSums(Z)
    Z <- log(Z)
    X <- log(X)

    y <- X %*% beta_star + rnorm(n, sd = sigma)

    # Return data
    return(list(X = X, Z = Z, y = y, Sig_B = Sig_B))
}

generate_multinom_data <- function(n, p, beta_star, sigma, rho, theta = NA, type = "multinom",overdispersion=5e+3 ) {
    # n: number of rows
    # p: number of columns
    # beta_star: true beta
    # sigma: noise level
    # tau: sd of measurement eror
    # rho: correlation between columns of W
    # theta: mean of W

    if (anyNA(theta)) theta <- c(rep(log(0.5 * p), 5), rep(0, p - 5))
    W <- AR(n, p, rho) + rep(theta, each = n)
    X <- exp(W) / rowSums(exp(W))

    ## Count matrix simulated from community compositions
    if (type == "multinom") {
        for (i in 1:nrow(W)) {
            nSeq <- rnbinom(1, mu = 10000, size = 25)
            W[i, ] <- rmultinom(1, nSeq, prob = X[i, ])[, 1]
        }
    } else if (type == "dirmult") {
        for (i in 1:nrow(W)) {
            Q <- rdirichlet(1, alpha = overdispersion * X[i, ])
            nSeq <- rnbinom(1, prob = 0.01, size = 300 / 0.99)
            W[i, ] <- rmultinom(1, nSeq, prob = Q)[, 1]
        }
    } else {
        stop("type must be either multinom or dirmult")
    }

    W <- W + 0.5 ## recommendated by 2023 BKA paper
    Z <- W / rowSums(W)

    Z <- log(Z)
    X <- log(X)
    y <- X %*% beta_star + rnorm(n, sd = sigma)

    return(list(X = X, Z = Z, y = y))
}

# Estimate the Sig_B through Monte-Carlo simulation

MC_varB <- function(n, p, beta_star, sigma, rho, theta = NA, N_MK = 10000 %/% p, type = "multinom",overdispersion=5e+3) {
    varB <- matrix(NA, N_MK, p)
    for (mk in 1:N_MK) {
        tmp_data <- generate_multinom_data(n, p, beta_star, sigma, rho, theta = theta, type = type,overdispersion=overdispersion)
        B <- tmp_data$Z - tmp_data$X
        for (i in 1:ncol(B)) {
            varB[mk, i] <- var(B[, i])
        }
    }
    varB <- colMeans(varB)
    varB[1:5] <- mean(varB[1:5])
    varB[6:p] <- mean(varB[6:p])
    Sig_B <- diag(varB)
    return(Sig_B)
}


# Function to generate the covariance matrix with geometric decay
AR_covariance_matrix <- function(n, rho) {
    # Initialize an empty covariance matrix
    Sigma <- matrix(0, nrow = n, ncol = n)

    # Calculate the covariance values based on geometric decay
    for (i in 1:n) {
        for (j in 1:n) {
            distance <- abs(i - j)
            Sigma[i, j] <- rho^distance
        }
    }

    return(Sigma)
}
sym_covariance_matrix <- function(n, rho) {
    Sigma <- matrix(0.5, nrow = n, ncol = n)
    Sigma <- Sigma + diag(0.5, nrow = n, ncol = n)
    return(Sigma)
}


# Function to normalize X
normalize <- function(X, scale = TRUE, center = TRUE) {
    n <- nrow(X)
    # centering
    if (center) X <- apply(X, 2, function(x) x - mean(x))
    # scale columns of z to sqrt(n)
    if (scale) X <- apply(X, 2, function(x) sqrt(n) * x / sqrt(sum(x^2)))
    return(X)
}

# function to get multivariate Gaussian variable
AR <- function(n, p, rho) {
    z <- matrix(rnorm(n * p), n, p)
    X <- matrix(NA, n, p)
    X[, 1] <- z[, 1]
    for (j in 2:(p)) {
        X[, j] <- rho * X[, j - 1] + sqrt(1 - rho^2) * z[, j]
    }
    return(X)
}


#--------------------------------------------
# Evaluate data
#--------------------------------------------

PE <- function(beta, X, beta_star) {
    n <- dim(X)[1]
    error <- (beta - beta_star) %*% t(X) %*% X %*% (beta - beta_star) / n
    as.numeric(error)
}
l_inf <- function(beta, beta_star) {
    return(max(abs(beta - beta_star)))
}
SE <- function(beta, beta_star) {
    return(sum((beta - beta_star)^2))
}
Falsility <- function(beta, beta_star) {
    FN <- sum(beta == 0 & beta_star != 0)
    FP <- sum(beta != 0 & beta_star == 0)
    FPR <- FP / sum(beta_star == 0)
    FNR <- FN / sum(beta_star != 0)
    return(list(FN = FN, FP = FP, FPR = FPR, FNR = FNR))
}

eval_results <- function(beta, X, beta_star) {
    Fa <- Falsility(beta, beta_star)
    list(PE = PE(beta, X, beta_star), SE = SE(beta, beta_star), FN = Fa$FN, FP = Fa$FP, l_inf = l_inf(beta, beta_star), FPR = Fa$FPR, FNR = Fa$FNR)
}
