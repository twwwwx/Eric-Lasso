source("R/eic.R")
source("real_DA/RDAdata.R")
source("simulation/generate_data.R")
set.seed(123456)

cov <- read.csv("real_DA/COV.csv", header = T)
N_sim <- 10
K <- 5

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

# print("perturb Z, 0.2~5")
# Z <- X
# for (i in 1:n) {
#     U <- runif(p, min = 0.05, max = 20) ## measurement error factors
#     ## OR U=runif(p,min=0.01,max=100)
#     Z[i, ] <- perturb(X[i, ], U)
# }

#######  Eric-Lasso analysis
## either with (y,X,Z) or (y1,X,Z)
regressionlist <- c("(y,X,Z)", "(y1,X,Z)")
for (m in 1:2) {
    # data <- get_data(X = X, y = y, y1 = y1, y_mode = m, min_perturb = 0.05, max_perturb = 20, K = K)
    print(regressionlist[m])
    if (m == 1) {
        y <- as.matrix(y) ## BMI
        print("Using (y,Z) to inference (y,X)")
    } else {
        y <- as.matrix(y1) ## cov-adjusted BMI
        print("Using (y1, Z) to inference (y1, X)")
    }
    
    # logX <- log(X)
    # logZ <- log(Z)
    # ind <- sample(1:nrow(X), K - nrow(X) %% K)
    # X.p <- rbind(logX, logX[ind, ])
    # Z.p <- rbind(logZ, logZ[ind, ])
    # y.p <- rbind(y, as.matrix(y[ind, ]))
    # n.p <- nrow(X.p)

    model_list <- list()
    model_list[["Eric"]] <- c(TRUE, TRUE)
    model_list[["CoDA"]] <- c(TRUE, FALSE)
    model_list[["CoCo"]] <- c(FALSE, TRUE)
    model_list[["Vani"]] <- c(FALSE, FALSE)
    results_df <- data.frame(lam=numeric(N_sim), MSE = numeric(N_sim),MSEout=numeric(N_sim), MAE=numeric(N_sim),MAEout=numeric(N_sim), selected=numeric(N_sim), sum_beta = numeric(N_sim))
    support <- list()
    write.table(t(c("Model",colnames(results_df),'p_value')),
        file = "RDAlog/RDAresults.csv", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
    )
    for (model_name in names(model_list)) {
        print(model_name)
        start_time <- Sys.time()
        set.seed(123456)
        betas <- matrix(NA, nrow = N_sim, ncol = p)
        intercepts <- matrix(NA, nrow = N_sim, ncol = 1)
        N_bs <- nrow(X) / 2
        bs_inds <- matrix(NA, nrow = N_sim, ncol = N_bs - N_bs %% K)
        for (i in 1:N_sim) {
            if (i %% 10 == 0) print(paste("Round", i))
            data <- get_data(X = X, y = y, y1 = y1, y_mode = m, min_perturb = 0.05, max_perturb = 20, K = K)
            # data <- get_norm_data(X = X, y = y, y1 = y1, y_mode = m, min_perturb = 0.1, max_perturb = 10, K = K)
            ## bootstrap
            bs_ind <- sample(1:N_bs, N_bs - N_bs %% K, replace = TRUE)
            bs_inds[i, ] <- bs_ind
            X. <- data$X[bs_ind, ]
            Z. <-data$Z[bs_ind, ]
            y. <- as.matrix(data$y[bs_ind, ])
            n. <- nrow(X.)

            cen_Z <- normalize(Z., scale = FALSE)
            cen_X <- normalize(X., scale = FALSE)
            B <- cen_Z - cen_X
            Sig_B <- t(B) %*% B / nrow(Z.)
            # Sig_B <- t(cen_Z) %*% cen_Z / nrow(Z.) - t(cen_X) %*% cen_X / nrow(X.)
            fit <- eic(Z = Z., y = y., n = n., p = p, scale.Z = FALSE, scale.y = FALSE, step = 100, K = K, mu = 10, earlyStopping_max = 30, Sig_B = Sig_B, etol = 1e-4, noise = "additive", penalty = "lasso", proj = model_list[[model_name]][2], constrain = model_list[[model_name]][1])
            betas[i, ] <- fit$beta.opt
            
            intercept <- fit$mean.y - fit$mean.Z %*% fit$beta.opt
            intercepts[i, ] <- intercept
            y_pred <- X. %*% fit$beta.opt + rep(intercept, each = n.)
            results_df$sum_beta[i] <- sum(fit$beta.opt)
            results_df$selected[i] <- sum(fit$beta.opt != 0)
            MSE <- mean((y_pred - y.)^2)
            MAE <- mean(abs(y_pred - y.))
            results_df$MSE[i] <- MSE
            results_df$MAE[i] <- MAE
            results_df$lam[i] <- fit$lambda.opt

        }
        print(Sys.time() - start_time)
        MSEout <- leave_one_out_MSE(X=data$X,y=data$y,sample_inds=bs_inds,betas=betas,intercepts=intercepts)
        MAEout <- leave_one_out_MAE(X=data$X,y=data$y,sample_inds=bs_inds,betas=betas,intercepts=intercepts)

        results_df$MSEout <- rep(MSEout,N_sim)
        results_df$MAEout <- rep(MAEout,N_sim)


        # save results
        support[[model_name]] <- betas
        p_value <- t.test(results_df$sum_beta, mu = 0)$p.value
        bootstrap_mean <- apply(results_df, 2, mean)
        bootstrap_mean_std <- apply(results_df, 2, sd) / sqrt(N_sim)
        evals <- rbind(bootstrap_mean, bootstrap_mean_std)
        print(model_name)
        print(evals)
        # save results
        l <- length(colnames(results_df))
        values <- sapply(c(1:ncol(evals))[-l], function(i) paste0(round(evals[1, i], 3), "(", round(evals[2, i], 3), ")"))
        sum_value <- paste0(bootstrap_mean[l], "(", bootstrap_mean_std[l], ")")
        values <- t(c(model_name,values, sum_value, p_value))
        write.table(values,
            file = "RDAlog/RDAresults.csv", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
        )
    }
    save(support, file = paste0("RDAlog/", "RDA", m, "_betas.RData"))
    cat("\n", file = "RDAlog/RDAresults.csv", append = TRUE)
}
