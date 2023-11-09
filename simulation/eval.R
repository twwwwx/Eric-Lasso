source("R/eic.R")
source("simulation/generate_data.R")
library("boot")


# --------------------------

# simulation settings
data_type <- "lognormal"
# data_type <- "dirichlet"
# data_type <- "multinom"
# data_type <- "dirmult"
N_sim <- 100
# n <- 50
# p <- 50
# create a list of different n and p values
np_list <- list(c(100, 100))

sigma <- 0.5
rho <- 0.5
tau <- 1.5


# model settings
model_list <- list()
model_list[["Eric"]] <- c(TRUE, TRUE)
# model_list[["CoDA"]] <- c(TRUE, FALSE)
# model_list[["CoCo"]] <- c(FALSE, TRUE)
# model_list[["Vani"]] <- c(FALSE, FALSE)
# --------------------------
file_name <- "results/results_table.csv"
file_name_sum <- "results/results_sum.csv"
colname <- t(c("model", "data_type", "n", "p", "N_sim", "tau", "rho", "lam", "SE", "PE", "l_inf", "FPR", "FNR"))
write.table(colname,
    file = file_name, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
)

colname <- t(c("model", "data_type", "n", "p", "N_sim", "tau", "rho", "sum", "p value"))
write.table(colname,
    file = file_name_sum, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
)
for (np in np_list) {
    n <- np[1]
    p <- np[2]
    beta_star <- c(1.2, -0.8, 0.7, 0, 0, -1.5, -1, 1.4, rep(0, p - 8))
    theta <- c(rep(log(0.2 * p), 5), rep(0, p - 5))
    for (i in seq_along(model_list)) {
        model <- model_list[[i]]
        model_name <- names(model_list)[i]

        constrain <- model[1]
        proj <- model[2]
        subdir <- "results/"
        if (data_type == "multinom" || data_type == "dirmult") {
            Sig_B_estimated <- MC_varB(n, p, beta_star, sigma, rho, theta = theta, N_MK = 1000000 %/% p, type = data_type, overdispersion = 5e+3)
            tau <- sqrt(Sig_B_estimated[6, 6])
            print(paste("estimated tau is", sqrt(Sig_B_estimated[6, 6])))
        }
        subdir_name <- paste0(subdir, model_name, "+", data_type, "_n", n, "_p", p, "_tau", tau, "_rho", rho, "_sigma", sigma, "_Nsim", N_sim)

        # --------------------------
        # run lasso code
        start_time <- Sys.time()
        set.seed(1234567)

        results_df <- data.frame(lambda = numeric(N_sim), SE = numeric(N_sim), PE = numeric(N_sim), linf = numeric(N_sim), FPR = numeric(N_sim), FNR = numeric(N_sim), sum_beta = numeric(N_sim))
        for (i in 1:N_sim) {
            if (i %% 20 == 0) {
                print(paste("Round", i))
            }
            # generate data

            if (data_type == "multinom" || data_type == "dirmult") {
                data <- generate_multinom_data(n, p, beta_star, sigma, rho, theta = theta, type = data_type, overdispersion = 5e+3)
                data$Sig_B <- Sig_B_estimated
            } else {
                data <- generate_data(n, p, beta_star, sigma, tau, rho, theta, type = data_type)
            }

            # record time
            fit_additive <- eic(Z = data$Z, y = data$y, n = n, p = p, scale.Z = FALSE, scale.y = FALSE, step = 100, K = 5, mu = 10, earlyStopping_max = 10, Sig_B = data$Sig_B, etol = 1e-4, noise = "additive", penalty = "lasso", proj = proj, constrain = constrain)

            # naive evaluation
            results_df$lambda[i] <- fit_additive$lambda.opt

            centered_X <- data$X - rep(fit_additive$mean.Z, each = n)
            measures <- eval_results(fit_additive$beta.opt, centered_X, beta_star)
            results_df$SE[i] <- measures$SE
            results_df$PE[i] <- measures$PE
            # results_df$FN[i] <- measures$FN
            # results_df$FP[i] <- measures$FP
            results_df$FPR[i] <- measures$FPR
            results_df$FNR[i] <- measures$FNR
            results_df$linf[i] <- measures$l_inf
            results_df$sum_beta[i] <- sum(fit_additive$beta.opt)
        }
        runtime <- Sys.time() - start_time

        #--------------------------
        # bootstrap and p value

        bootstrap_mean <- apply(results_df, 2, mean)
        bootstrap_mean_std <- apply(results_df, 2, sd) / sqrt(N_sim)
        p_value <- t.test(results_df$sum_beta, mu = 0)$p.value

        if (model_name == "coco" && data_type == "dirichlet") {
            tmp <- results_df$sum_beta
            save(tmp, file = paste0(subdir_name, ".RData"))
        }


        #--------------------------
        # save results


        lam_value <- round(bootstrap_mean[1], 2)
        SE_value <- paste0(round(bootstrap_mean[2], 3), "(", round(bootstrap_mean_std[2], 3), ")")
        PE_value <- paste0(round(bootstrap_mean[3], 3), "(", round(bootstrap_mean_std[3], 3), ")")
        l_inf_value <- paste0(round(bootstrap_mean[4], 3), "(", round(bootstrap_mean_std[4], 3), ")")
        FPR_value <- paste0(round(bootstrap_mean[5], 3), "(", round(bootstrap_mean_std[5], 3), ")")
        FNR_value <- paste0(round(bootstrap_mean[6], 3), "(", round(bootstrap_mean_std[6], 3), ")")
        values <- t(as.matrix(c(model_name, data_type, n, p, N_sim, tau, rho, lam_value, SE_value, PE_value, l_inf_value, FPR_value, FNR_value)))
        sum_value <- paste0(round(bootstrap_mean[7], 3), "(", round(bootstrap_mean_std[7], 3), ")")
        sum_values <- t(as.matrix(c(model_name, data_type, n, p, N_sim, tau, rho, sum_value, p_value)))

        write.table(values,
            file = file_name, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
        )
        write.table(sum_values,
            file = file_name_sum, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
        )

        # print results
        evals <- rbind(bootstrap_mean, bootstrap_mean_std)
        colnames(evals) <- colnames(results_df)
        print(subdir_name)
        print(evals)

        print(runtime)
    }
    cat("\n", file = file_name, append = TRUE)
    cat("\n", file = file_name_sum, append = TRUE)
}
