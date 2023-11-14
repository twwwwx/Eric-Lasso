
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

MSE <- function(beta, X, y) {
    error <- mean((y - X %*% beta)^2)
    as.numeric(error)
}
MAE <- function(beta, X, y) {
    error <- mean(abs(y - X %*% beta))
    as.numeric(error)
}