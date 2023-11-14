source("R/eic.R")
source("real_DA/RDAdata.R")
cov <- read.csv("real_DA/COV.csv", header = T)
N_sim <- 100
K <- 5
set.seed(123456)


ltn <- function(X) {
    return(matrix(as.numeric(unlist(X)), nrow = nrow(X), ncol = ncol(X)))
}
cov <- ltn(cov)
y <- cov[, 1] ## BMI
fit <- lm(y ~ cov[, 2:3])
y1 <- fit$residuals ## cov-adjusted BMI

otu <- read.csv("real_DA/OTU.csv", header = T)
otu <- otu + 0.5 ## followed suggestion of BKA2022 paper
X <- otu
for (i in 1:nrow(X)) {
    X[i, ] <- X[i, ] / sum(otu[i, ])
}
# apply(X, 1, sum)
X <- ltn(X)
n <- nrow(X)
p <- ncol(X)
approx <- numeric(10)
no_approx <- numeric(10)
MSE <- numeric(10)

proj <- FALSE
constrain <- TRUE

for(i in 1:100){
data <- get_norm_data(X,y,y1,1,0.05,20)
# data <- get_data(X,y,y1,1,1/4,4)
Z. <- data$Z[-1,]
X. <- data$X[-1,]
X.norm <- normalize(X.,cen=T,scale=F)
Z.norm <- normalize(Z.,cen=T,scale=F)
B. <- Z.norm - X.norm
y. <- as.matrix(data$y[-1,])
n. <- nrow(Z.)
Sig_B <- t(B.) %*% B. / nrow(B.)
print(range(diag(Sig_B)))
# fit <- eic(Z = Z., y = y., n = n., p = p, scale.Z = FALSE, scale.y = FALSE, step = 100, K = 5, mu = 10, earlyStopping_max = 30, Sig_B = Sig_B, etol = 1e-4, noise = "additive", penalty = "lasso", proj = proj, constrain = constrain)
delta <- t(Z.norm) %*% Z.norm / nrow(Z.) - t(X.norm) %*% X. / nrow(X.norm)
est1 <- t(B.) %*% B. / nrow(B.)

# intercept <- fit$mean.y - fit$mean.Z %*% fit$beta.opt
# y_pred <- X. %*% fit$beta.opt + rep(intercept, each = n.)
# MSE[i] <- mean((y_pred - y.)^2)
approx[i] <- mean(abs(delta - est1))
no_approx[i] <- mean(abs(delta))
# break
}
print(mean(approx))
print(mean(no_approx))
# print(mean(MSE))
