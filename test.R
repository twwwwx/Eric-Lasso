source("simulation/generate_data.R")
p <- 8
beta_star <- c(1.2, -0.8, 0.7, 0, 0, -1.5, -1, 1.4, rep(0, p - 8))
theta <- c(rep(log(0.5 * p), 5), rep(0, p - 5))
n <- 10
sigma <- 0.5
rho <- 0.5

data <- generate_multinom_data(n, p, beta_star, sigma, rho, theta = theta)
sig <- MC_varB(n, p, beta_star, sigma, rho, theta = theta, N_MK=10000 %/% p)
data
