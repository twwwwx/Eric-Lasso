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
    }
    Z <- X * exp(B)
    Z <- Z / rowSums(Z)
    Z <- log(Z)
    X <- log(X)

    y <- X %*% beta_star + rnorm(n, sd = sigma)

    # Return data
    return(list(X = X, Z = Z, y = y, Sig_B = Sig_B))
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
