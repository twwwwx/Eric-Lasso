source("R/eic.R")
source("simulation/generate_data.R")
N_sim <- 4

##### read data #####
ymat <- read.csv("real_DA/responses303.csv", header = T)
otu1 <- read.csv("real_DA/pouchmat303.csv", header = T)
otu2 <- read.csv("real_DA/PPImat303.csv", header = T)

y1 <- ymat[, 1]
y2 <- ymat[, 2]

delta1 <- otu1 != 0
cutoff <- min(otu1[otu1 > 0], otu2[otu2 > 0])

for (i in 1:nrow(otu1)) {
    for (j in 1:ncol(otu1)) {
        if (otu1[i, j] == 0) otu1[i, j] <- cutoff
    }
}
for (i in 1:nrow(otu1)) {
    libsize <- sum(otu1[i, ])
    otu1[i, ] <- otu1[i, ] / libsize
}


delta2 <- otu2 != 0
for (i in 1:nrow(otu2)) {
    for (j in 1:ncol(otu2)) {
        if (otu2[i, j] == 0) otu2[i, j] <- cutoff
    }
}
for (i in 1:nrow(otu2)) {
    libsize <- sum(otu2[i, ])
    otu2[i, ] <- otu2[i, ] / libsize
}


#######  Eric-Lasso analysis
#####  Using (y1,otu2) to inference (y1,otu1)
## or using (y2, otu1) to inference (y2, otu2)
regressionlist <- c("y1~otu2", "y2~otu1")
for(i in 1:2){
    print(regressionlist[i])
    if(i==1){
        print("Using (y1,otu2) to inference (y1,otu1)")
        y <- as.matrix(y1)
        X <- as.matrix(log(otu1))
        Z <- as.matrix(log(otu2))
    }else{
        print("Using (y2, otu1) to inference (y2, otu2)")
        y <- as.matrix(y2)
        X <- as.matrix(log(otu2))
        Z <- as.matrix(log(otu1))
    }

    print("apply scale to X")
    X <- apply(X, 2, function(x) (x-mean(x))/sd(x))
    p <- ncol(X)
    print("apply scale to Z")
    Z <- apply(Z, 2, function(x) (x-mean(x))/sd(x))
    B <- Z - X
    print("compute Sig_B by var(B)")
    Sig_B <- t(B) %*% B / nrow(B)



    # # pad the matrix for 5-fold CV
    set.seed(123456)
    ind <- sample(1:nrow(X), 5 - nrow(X) %% 5)
    X.p <- rbind(X, X[ind, ])
    Z.p <- rbind(Z, Z[ind, ])
    y.p <- rbind(y, as.matrix(y[ind, ]))
    n.p <- nrow(X.p)

    # ## get the true beta using coda regression
    # fit_t <- eic(Z = X.p, y = y.p, n = n.p, p = p, scale.Z = FALSE, scale.y = FALSE, step = 100, K = 5, mu = 10, earlyStopping_max = 10, Sig_B = NULL, etol = 1e-4, noise = "additive", penalty = "lasso", proj = FALSE, constrain = TRUE)
    # beta_star <- fit_t$beta.opt

    model_list <- list()
    model_list[["Eric"]] <- c(TRUE, TRUE)
    model_list[["CoDA"]] <- c(TRUE, FALSE)
    model_list[["CoCo"]] <- c(FALSE, TRUE)
    model_list[["Vani"]] <- c(FALSE, FALSE)
    i=1
    results_df <- data.frame(model = character(N_sim), lambda = numeric(N_sim), SE = numeric(N_sim), PE = numeric(N_sim), linf = numeric(N_sim), FPR = numeric(N_sim), FNR = numeric(N_sim), sum_beta = numeric(N_sim),MSE=numeric(N_sim))
    for (model_name in names(model_list)) {
        start_time <- Sys.time()
        set.seed(123456)
        fit <- eic(Z = Z.p, y = y.p, n = n.p, p = p, scale.Z = FALSE, scale.y = FALSE, step = 100, K = 5, mu = 10, earlyStopping_max = 30, Sig_B = Sig_B, etol = 1e-4, noise = "additive", penalty = "lasso", proj = model_list[[model_name]][2], constrain = model_list[[model_name]][1])
        # centered_X <- X.p - rep(fit$mean.Z, each = n.p)
        # measures <- eval_results(fit$beta.opt, centered_X, beta_star)
        results_df$model[i] <- model_name
        # results_df$SE[i] <- measures$SE
        # results_df$PE[i] <- measures$PE
        results_df$FPR[i] <- sum(fit$beta.opt != 0)
        # results_df$FNR[i] <- measures$FNR
        # results_df$linf[i] <- measures$l_inf

        intercept <- fit$mean.y - fit$mean.Z %*% fit$beta.opt
        y_pred <- X.p %*% fit$beta.opt + rep(intercept, each = n.p)
        results_df$sum_beta[i] <- sum(fit$beta.opt)
        results_df$lambda[i] <- fit$lambda.opt
        results_df$MSE[i] <- mean((y_pred - y.p)^2)
        print(results_df[i,])  
        print(Sys.time() - start_time)
        i=i+1
    }
    print(results_df)  
}