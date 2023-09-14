source("R/eic.R")
source("simulation/generate_data.R")
library("boot")


# --------------------------

# simulation settings

N_sim <- 50
N_bs <- NA
n <- 100
p <- 200
sigma <- 0.5
rho <- 0.2
id <- FALSE
normalize <- TRUE
model_list <- list()
model_list[["eic"]] <- c(TRUE, TRUE)
model_list[["coda"]] <- c(TRUE, FALSE)
model_list[["coco"]] <- c(FALSE, TRUE)
# --------------------------
for (tau in c(0.25, 0.5, 0.75)) {
    for (i in seq_along(model_list)) {
        model <- model_list[[i]]
        model_name <- names(model_list)[i]
        if (normalize) {
            model_name <- paste0(model_name, "+norm")
        }
        if (id) {
            model_name <- paste0(model_name, "+id")
        }
        constrain <- model[1]
        proj <- model[2]
        subdir <- paste0("results/", model_name, "/")
        subdir_name <- paste0(subdir, "_n", n, "_p", p, "_tau", tau, "_rho", rho, "_sigma", sigma, "_id", id, "_Nsim", N_sim, "_Nbs", N_bs, "_norm", normalize)

        # --------------------------
        # run lasso code
        start_time <- Sys.time()
        set.seed(1234567)

        results_df <- data.frame(lambda = numeric(N_sim), SE = numeric(N_sim), PE = numeric(N_sim), FN = numeric(N_sim), FP = numeric(N_sim), linf = numeric(N_sim), FPR = numeric(N_sim), FNR = numeric(N_sim))
        for (i in 1:N_sim) {
            if (i %% 20 == 0) {
                print(paste("Round", i))
            }
            # generate data
            beta_star <- c(1.2, -0.8, 0.7, 0, 0, -1.5, -1, 1.4, rep(0, p - 8))
            theta <- c(rep(log(0.5 * p), 5), rep(0, p - 5))
            data <- generate_data(n, p, beta_star, sigma, tau, rho, theta, id = id, normalize = normalize)

            # record time
            fit_additive <- eic(Z = data$Z, y = data$y, n = n, p = p, scale.Z = FALSE, scale.y = FALSE, step = 100, K = 5, mu = 10, earlyStopping_max = 10, Sig_B = data$Sig_B, etol = 1e-4, noise = "additive", penalty = "lasso", proj = proj, constrain = constrain)

            # naive evaluation
            results_df$lambda[i] <- fit_additive$lambda.opt

            centered_X <- data$X - rep(fit_additive$mean.Z, each = n)
            measures <- eval_results(fit_additive$beta, centered_X, beta_star)
            results_df$SE[i] <- measures$SE
            results_df$PE[i] <- measures$PE
            results_df$FN[i] <- measures$FN
            results_df$FP[i] <- measures$FP
            results_df$FPR[i] <- measures$FPR
            results_df$FNR[i] <- measures$FNR
            results_df$linf[i] <- measures$l_inf
        }
        runtime <- Sys.time() - start_time

        #--------------------------
        # bootstrap

        if(!is.na(N_bs)){
            med_stats <- function(data, i) {
                apply(data[i, ], 2, median)
            }
            bootstrap_results <- boot(data = results_df, statistic = med_stats, R = N_bs)

            bootstrap_median <- apply(bootstrap_results$t, 2, mean)
            bootstrap_median_std <- apply(bootstrap_results$t, 2, sd)
        }
        else{
            bootstrap_median <- apply(results_df, 2, mean)
            bootstrap_median_std <- apply(results_df, 2, sd)
        }

        #--------------------------
        # save results
        file_name <- "results/results_table.csv"
        if (!file.exists(file_name)) {
            write.table(c("model", "n", "p", "N_sim", "N_bs", "tau", "rho", "lam", "SE", "PE", "FN", "FP", "l_inf", "FPR", "FNR"),
                file = file_name, sep = ",", row.names = FALSE, col.names = FALSE
            )
        }
        lam_value <- round(bootstrap_median[1], 2)
        SE_value <- paste0(round(bootstrap_median[2], 2), "(", round(bootstrap_median_std[2], 2), ")")
        PE_value <- paste0(round(bootstrap_median[3], 2), "(", round(bootstrap_median_std[3], 2), ")")
        FN_value <- paste0(round(bootstrap_median[4], 2), "(", round(bootstrap_median_std[4], 2), ")")
        FP_value <- paste0(round(bootstrap_median[5], 2), "(", round(bootstrap_median_std[5], 2), ")")
        l_inf_value <- paste0(round(bootstrap_median[6], 2), "(", round(bootstrap_median_std[6], 2), ")")
        FPR_value <- paste0(round(bootstrap_median[7], 2), "(", round(bootstrap_median_std[7], 2), ")")
        FNR_value <- paste0(round(bootstrap_median[8], 2), "(", round(bootstrap_median_std[8], 2), ")")
        values <- t(as.matrix(c(model_name, n, p, N_sim, N_bs, tau, rho, lam_value, SE_value, PE_value, FN_value, FP_value, l_inf_value, FPR_value, FNR_value)))

        write.table(values,
            file = file_name, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
        )

        # print results
        evals <- rbind(bootstrap_median, bootstrap_median_std)
        colnames(evals) <- colnames(results_df)
        print(subdir_name)
        print(evals)

        print(runtime)
    }
}
