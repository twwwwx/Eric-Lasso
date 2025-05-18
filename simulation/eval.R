source("R/eic.R")
source("simulation/generate_data.R")
source("simulation/functions.R")
library("boot")


# --------------------------

# simulation settings
# data_type <- "lognormal"
data_type <- "dirichlet"
# data_type <- "multinom"
# data_type <- "dirmult"
# data_type_list <- c("dirichlet", "dirmult")
# data_type_list <- c("dirmult")
N_sim <- 100
# create a list of different n and p values
np_list <- list(c(100,200))
# np_list <- list(c(100,200),c(250,400),c(500,500),c(550,700))

sigma <- 0.5
rho <- 0.5
tau <- 0.5
tau_list <- seq(0.1, 1.9, 0.2)
overdispersion = 5e+3
# model settings
model_list <- list()
# model_list[["Debi"]] <- c(FALSE)
model_list[["Eric"]] <- c(TRUE, TRUE)
# model_list[["Coda"]] <- c(TRUE, FALSE)
# model_list[["CoCo"]] <- c(FALSE, TRUE)
# model_list[["Vani"]] <- c(FALSE, FALSE)
# -------------------------- 
file_name <- "results/results_table.csv"
file_name_sum <- "results/results_sum.csv"
colname <- t(c("model", "data_type", "n", "p", "N_sim", "tau", "rho", "lam", "SE", "PE", "l_inf", "FPR", "FNR","TPR"))
write.table(colname,
    file = file_name, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
)

colname <- t(c("model", "data_type", "n", "p", "N_sim", "tau", "rho", "sum", "p value"))
write.table(colname,
    file = file_name_sum, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
)
for (np in np_list) {
# for (tau in tau_list) {
# for (data_type in data_type_list) {
    # np = c(100,100)
    n <- np[1]
    p <- np[2]
    beta_star <- c(1.2, -0.8, 0.7, 0, 0, -1.5, -1, 1.4, rep(0, p - 8))
    theta <- c(rep(log(0.2 * p), 5), rep(0, p - 5))
    # theta = rep(0, p)
    for (i in seq_along(model_list)) {
        model <- model_list[[i]]
        model_name <- names(model_list)[i]
        if(length(model)== 2){
            constrain <- model[1]
            proj <- model[2]
        }
        subdir <- "data/"
        if (data_type == "multinom" || data_type == "dirmult") {

            ### !!! modify the overdispersion
            Sig_B_estimated <- MC_varB(n, p, beta_star, sigma, rho, theta = theta, N_MK = 100000 %/% p, type = data_type, overdispersion = overdispersion)
            # Sig_B_estimated <- MC_varB(n, p, beta_star, sigma, rho, theta = theta, N_MK = 100000 %/% p, type = data_type, overdispersion = 5e+3)
            tau <- sqrt(Sig_B_estimated[6, 6])
            print(paste("estimated tau is", sqrt(Sig_B_estimated[6, 6])))
        }
        subdir_name <- paste0(subdir, model_name, "+", data_type, "_n", n, "_p", p, "_tau", tau, "_rho", rho, "_sigma", sigma, "_Nsim", N_sim)

        # --------------------------
        # run lasso code
        start_time <- Sys.time()
        set.seed(1234567)

        results_df <- data.frame(lambda = numeric(N_sim), SE = numeric(N_sim), PE = numeric(N_sim), linf = numeric(N_sim), FPR = numeric(N_sim),FNR = numeric(N_sim), TPR = numeric(N_sim), sum_beta = numeric(N_sim))
        results_bias <- matrix(NA, nrow = N_sim, ncol = p)
        for (i in 1:N_sim) {
            # i = i + 20
            if (i %% 2 == 0) {
                print(paste("Round", i))
            }
            # generate data
            if (data_type == "multinom" || data_type == "dirmult") {
                data <- generate_multinom_data(n, p, beta_star, sigma, rho, theta = theta, type = data_type, overdispersion = overdispersion)
                data$Sig_B <- Sig_B_estimated
            } else {
                data <- generate_data(n, p, beta_star, sigma, tau, rho, theta, type = data_type)
            }

            # record time
            if(length(model)== 2){
                # Note: if p > 700, step is set to be 50
                fit_additive <- eic(Z = data$Z, y = data$y, n = n, p = p, scale.Z = FALSE, scale.y = FALSE, step = 100, K = 5, mu = 10, earlyStopping_max = 10, Sig_B = data$Sig_B, etol = 1e-4, 
                                    noise = "additive",
                                    penalty = "lasso",
                                    proj = proj, 
                                    constrain = constrain)
            }else {
                Zp <- ref_transform(data$Z)
                mu_u = rep(0,p)
                Sigma_x_lat = AR_covariance_matrix(p, rho)
                Sigma_w_lat = Sigma_x_lat+data$Sig_B
                supp_data = list(mu_x = theta, Sigma_x = Sigma_x_lat, Sigma_w = Sigma_w_lat, mu_w = theta)
                # fit_additive <- fn_proposed(VV = Zp,y = data$y, alpha_real=beta_star, W = exp(data$Z),mu_u = mu_u, Sigma_u = data$Sig_B,EstimateSigma = T, Noestimate.data = supp_data)
                fit_additive <- fn_proposed(VV = Zp,y = data$y, alpha_real=beta_star, W = exp(data$Z),mu_u = mu_u, Sigma_u = data$Sig_B,EstimateSigma = T, Noestimate.data = supp_data)
                
                # fit_additive$beta.opt <- fit_additive$beta.test
                # print(fit_additive$beta.opt[1:10])
                # print(beta_star[1:10])
            }   
            # naive evaluation
            results_df$lambda[i] <- fit_additive$lambda.opt

            centered_X <- data$X - rep(fit_additive$mean.Z, each = n)
            # centered_X_ref <- ref_transform(centered_X)
            # measures <- eval_results(fit_additive$beta.opt[1:(p-1)], centered_X_ref, beta_star[1:(p-1)], beta_test = fit_additive$beta.test[1:(p-1)])
            # measures <- eval_results(fit_additive$beta.lasso, centered_X, beta_star)
            measures <- eval_results(fit_additive$beta.opt, centered_X, beta_star, beta_test = fit_additive$beta.test)
            results_df$SE[i] <- measures$SE
            results_df$PE[i] <- measures$PE
            results_df$FPR[i] <- measures$FPR
            results_df$FNR[i] <- measures$FNR
            results_df$TPR[i] <- 1 - measures$FNR
            results_df$linf[i] <- measures$l_inf
            results_df$sum_beta[i] <- sum(fit_additive$beta.opt)

            results_bias[i,] <- fit_additive$beta.opt - beta_star
        }
        runtime <- Sys.time() - start_time

        #--------------------------
        # bootstrap and p value
        
        # bootstrap_metric <- function(data, R = 500, N = nrow(data)) {
            # median_fn <- function(data, indices) median(data[indices])
        #     
        #     boot_results <- boot(data = data, statistic = median_fn, R = R)
        #     
        #     list(
        #         bootstrap_median = mean(boot_results$t),
        #         bootstrap_se = sd(boot_results$t)
        #     )
        # }

        # Apply bootstrap function to each column
        # bootstrap_results <- lapply(results_df, bootstrap_metric)

        # bootstrap_mean = sapply(bootstrap_results, function(x) x$bootstrap_median) 
        # bootstrap_mean_std = sapply(bootstrap_results, function(x) x$bootstrap_se)
        bootstrap_mean <- apply(results_df, 2, mean)
        # average_bias <- apply(results_bias, 2, mean)
        bootstrap_mean_std <- apply(results_df, 2, sd) / sqrt(N_sim)
        # average_bias_std <- apply(results_bias, 2, sd) / sqrt(N_sim)
        p_value <- t.test(results_df$sum_beta, mu = 0)$p.value
        # print(results_bias[c(15:20),])
        # print(average_bias)
        # print(average_bias_std)
        # if (model_name == "CoCo" && data_type == "dirichlet") {
        #     tmp <- results_df$sum_beta
        #     save(tmp, file = paste0(subdir_name, ".RData"))
        # }


        #--------------------------
        # save results


        lam_value <- round(bootstrap_mean[1], 2)
        SE_value <- paste0(round(bootstrap_mean[2], 3), "(", round(bootstrap_mean_std[2], 3), ")")
        PE_value <- paste0(round(bootstrap_mean[3], 3), "(", round(bootstrap_mean_std[3], 3), ")")
        l_inf_value <- paste0(round(bootstrap_mean[4], 3), "(", round(bootstrap_mean_std[4], 3), ")")
        FPR_value <- paste0(round(bootstrap_mean[5], 3), "(", round(bootstrap_mean_std[5], 3), ")")
        FNR_value <- paste0(round(bootstrap_mean[6], 3), "(", round(bootstrap_mean_std[6], 3), ")")
        TPR_value <- paste0(round(bootstrap_mean[7], 3), "(", round(bootstrap_mean_std[7], 3), ")")
        values <- t(as.matrix(c(model_name, data_type, n, p, N_sim, tau, rho, lam_value, SE_value, PE_value, l_inf_value, FPR_value, FNR_value,TPR_value)))
        sum_value <- bootstrap_mean[8]
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
        # print(results_df$SE)
        # print(results_df$PE)
    }
    cat("\n", file = file_name, append = TRUE)
    cat("\n", file = file_name_sum, append = TRUE)
}
# save(results_bias, file = "bias_100_200_Debi.RData")
