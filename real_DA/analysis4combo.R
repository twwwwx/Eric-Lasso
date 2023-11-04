source("simulation/generate_data.R")
source("R/eic.R")

cov <- read.csv("real_DA/COV.csv", header = T)
N_sim <- 100
K <- 5
dim(cov)
ltn <- function(X) {
    return(matrix(as.numeric(unlist(X)), nrow = nrow(X), ncol = ncol(X)))
}
cov <- ltn(cov)
y <- cov[, 1] ## BMI
fit <- lm(y ~ cov[, 2:3])
y1 <- fit$residuals ## cov-adjusted BMI
length(y)

otu <- read.csv("real_DA/OTU.csv", header = T)
dim(otu) ## n=96, p=80
range(otu)
otu <- otu + 0.5 ## followed suggestion of BKA2022 paper
X <- otu
for (i in 1:nrow(X)) {
    X[i, ] <- X[i, ] / sum(otu[i, ])
}
# apply(X, 1, sum)
X <- ltn(X)

perturb <- function(a, b) {
    p <- length(a)
    out <- rep(NA, p)
    for (i in 1:p) {
        out[i] <- a[i] * b[i] / t(a) %*% b
    }
    return(out)
}

Z <- X
n <- nrow(Z)
p <- ncol(Z)
for (i in 1:n) {
    U <- runif(p, min = 0.1, max = 10) ## measurement error factors
    ## OR U=runif(p,min=0.01,max=100)
    Z[i, ] <- perturb(X[i, ], U)
}
## From here you can perform Eric-Lasso analysis



#######  Eric-Lasso analysis
## either with (y,X,Z) or (y1,X,Z)
regressionlist <- c("(y,X,Z)", "(y1,X,Z)")
for (m in 1:2) {
    print(regressionlist[m])
    if (m == 1) {
        print("Using (y,Z) to inference (y,X)")
        y <- as.matrix(y)
        X <- as.matrix(X)
        Z <- as.matrix(Z)
    } else {
        print("Using (y1, Z) to inference (y1, X)")
        y <- as.matrix(y1)
        X <- as.matrix(X)
        Z <- as.matrix(Z)
    }

    # print("apply scale to X")
    # logX <- apply(log(X), 2, function(x) (x - mean(x)) / sd(x))
    logX <- log(X)
    # print("apply scale to Z")
    # logZ <- apply(log(Z), 2, function(x) (x - mean(x)) / sd(x))
    logZ <- log(Z)
    set.seed(123456)
    ind <- sample(1:nrow(X), K - nrow(X) %% K)
    X.p <- rbind(logX, logX[ind, ])
    Z.p <- rbind(logZ, logZ[ind, ])
    y.p <- rbind(y, as.matrix(y[ind, ]))
    n.p <- nrow(X.p)


    B <- Z.p - X.p
    print("compute Sig_B by var(B)")
    Sig_B <- t(B) %*% B / nrow(B)
    print("range of Sig_B")
    range(Sig_B)


    delta <- cov(logZ) - cov(logX)
    est1 <- Sig_B
    est2 <- -Sig_B
    print("Use Sig_B for approximation")
    print(mean(abs(delta - est1)))
    print("No approximation")
    print(mean(abs(delta)))

    # # pad the matrix for 5-fold CV

    beta_star <- eic(Z = X.p, y = y.p, n = n.p, p = p, scale.Z = FALSE, scale.y = FALSE, step = 100, K = K, mu = 10, earlyStopping_max = 10, Sig_B = NULL, etol = 1e-4, noise = "additive", penalty = "lasso", proj = FALSE, constrain = TRUE)$beta.opt


    # ## get the true beta using coda regression
    # fit_t <- eic(Z = X.p, y = y.p, n = n.p, p = p, scale.Z = FALSE, scale.y = FALSE, step = 100, K = 5, mu = 10, earlyStopping_max = 10, Sig_B = NULL, etol = 1e-4, noise = "additive", penalty = "lasso", proj = FALSE, constrain = TRUE)
    # beta_star <- fit_t$beta.opt

    model_list <- list()
    model_list[["Eric"]] <- c(TRUE, TRUE)
    model_list[["CoDA"]] <- c(TRUE, FALSE)
    model_list[["CoCo"]] <- c(FALSE, TRUE)
    model_list[["Vani"]] <- c(FALSE, FALSE)
    results_df <- data.frame(lambda = numeric(N_sim), SE = numeric(N_sim), PE = numeric(N_sim), linf = numeric(N_sim), FPR = numeric(N_sim), FNR = numeric(N_sim), sum_beta = numeric(N_sim), MSE = numeric(N_sim))
    for (model_name in names(model_list)) {
        print(model_name)
        start_time <- Sys.time()
        set.seed(123456)
        for (i in 1:N_sim) {
            if (i %% 10 == 0) print(paste("Round", i))
            ## bootstrapd
            bs_ind <- sample(1:nrow(X), nrow(X) - nrow(X) %% K, replace = TRUE)
            X. <- logX[bs_ind, ]
            Z. <- logZ[bs_ind, ]
            y. <- as.matrix(y[bs_ind, ])
            n. <- nrow(X.)
            fit <- eic(Z = Z., y = y., n = n., p = p, scale.Z = FALSE, scale.y = FALSE, step = 100, K = K, mu = K, earlyStopping_max = 30, Sig_B = Sig_B, etol = 1e-4, noise = "additive", penalty = "lasso", proj = model_list[[model_name]][2], constrain = model_list[[model_name]][1])
            centered_X <- X.p - rep(fit$mean.Z, each = n.p)
            measures <- eval_results(fit$beta.opt, centered_X, beta_star)
            # results_df$model[i] <- model_name
            results_df$lambda[i] <- fit$lambda.opt
            results_df$SE[i] <- measures$SE
            results_df$PE[i] <- measures$PE
            results_df$FPR[i] <- measures$FPR
            results_df$FNR[i] <- measures$FNR
            results_df$linf[i] <- measures$l_inf

            intercept <- fit$mean.y - fit$mean.Z %*% fit$beta.opt
            y_pred <- X.p %*% fit$beta.opt + rep(intercept, each = n.p)
            results_df$sum_beta[i] <- sum(fit$beta.opt)
            results_df$seleted[i] <- sum(fit$beta.opt != 0)
            results_df$MSE[i] <- mean((y_pred - y.p)^2)
        }
        print(Sys.time() - start_time)
        bootstrap_mean <- apply(results_df, 2, mean)
        bootstrap_mean_std <- apply(results_df, 2, sd) / sqrt(N_sim)
        evals <- rbind(bootstrap_mean, bootstrap_mean_std)
        print(model_name)
        print(evals)
    }
}
