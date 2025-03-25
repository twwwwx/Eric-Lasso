source("R/eic.R")
source("simulation/generate_data.R")
source("simulation/functions.R")
library("boot")


# --------------------------

# simulation settings
# data_type <- "lognormal"
# data_type <- "dirichlet"
# data_type <- "multinom"
# data_type <- "dirmult"
file_name <- "data/ROC_S2.csv"
N_sim <- 100
# create a list of different n and p values
# np_list <- list(c(100,100))
# fpr_list <- seq(0.4, 0.7, 0.05)    
# fpr_list <- c(seq(0, 0.002, 0.0005), seq(0.6, 0.8, 0.05))

sigma <- 0.5
rho <- 0.5
tau <- 2.5
# tau_list <- c(10.5,11.5,12.5)

# model settings
model_list <- list()
# model_list[["Debi"]] <- c(FALSE)
# model_list[["LassoII"]] <- c(TRUE)

model_list[["Eric"]] <- c(TRUE, TRUE)
model_list[["Coda"]] <- c(TRUE, FALSE)
# model_list[["CoCo"]] <- c(FALSE, TRUE)
# model_list[["Vani"]] <- c(FALSE, FALSE)
n <- 100
p <- 100


# -------------------------- 
file_name_sum <- "results/results_sum.csv"
colname <- t(c("model", "data_type", "n", "p", "N_sim", "tau", "rho", "lam", "SE", "PE", "l_inf", "FPR", "FNR","TPR"))
# write.table(colname,
#     file = file_name, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
# )

colname <- t(c("model", "data_type", "n", "p", "N_sim", "tau", "rho", "sum", "p value"))
# write.table(colname,
#     file = file_name_sum, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
# )
# for (fpr in fpr_list) {
beta_star <- c(1.2, -0.8, 0.7, 0, 0, -1.5, -1, 1.4, rep(0, p - 8))
theta <- c(rep(log(0.2 * p), 5), rep(0, p - 5))
data <- generate_data(n, p, beta_star, sigma, tau, rho, theta, type = data_type)
lambda_max <- lambda_max(Z = data$Z, y = data$y, n = n, ratio_matrix = NULL, noise = "additive")
lambda_min <- 0.01 * lambda_max
lambda_list <- emdbook::lseq(lambda_max, lambda_min, 20)
print(lambda_list)
# ----------------------------
# run simulation
for( lambda_step in lambda_list){
    for (i in seq_along(model_list)) {
        model <- model_list[[i]]
        model_name <- names(model_list)[i]
        if(length(model)== 2){
            constrain <- model[1]
            proj <- model[2]
        }
        subdir <- "data/"
        if (data_type == "multinom" || data_type == "dirmult") {
            Sig_B_estimated <- MC_varB(n, p, beta_star, sigma, rho, theta = theta, N_MK = 100000 %/% p, type = data_type, overdispersion = 5e+3)
            tau <- sqrt(Sig_B_estimated[6, 6])
            print(paste("estimated tau is", sqrt(Sig_B_estimated[6, 6])))
        }
        subdir_name <- paste0(subdir, model_name, "+", data_type, "_n", n, "_p", p, "_tau", tau, "_rho", rho, "_sigma", sigma, "_Nsim", N_sim)

        # --------------------------
        # run lasso code
        start_time <- Sys.time()
        set.seed(1234567)

        results_df <- data.frame(lambda = numeric(N_sim), SE = numeric(N_sim), PE = numeric(N_sim), linf = numeric(N_sim), FPR = numeric(N_sim),FNR = numeric(N_sim), TPR = numeric(N_sim), sum_beta = numeric(N_sim))
        results_bias <- matrix(NA, nrow = N_sim, ncol = 11)
        for (i in 1:N_sim) {
            if (i %% 10 == 0) {
                print(paste("Round", i))
            }
            # generate data

            if (data_type == "multinom" || data_type == "dirmult") {
                data <- generate_multinom_data(n, p, beta_star, sigma, rho, theta = theta, type = data_type, overdispersion = 5e+3)
                data$Sig_B <- Sig_B_estimated
            } else {
                data <- generate_data(n, p, beta_star, sigma, tau, rho, theta, type = data_type)
                print("data is generated")
            }

            # record time
            if(length(model)== 2){
                fit_additive <- eic_one_step(Z = data$Z, y = data$y, n = n, p = p, lambda=lambda_step, scale.Z = FALSE, scale.y = FALSE, step = 50, K = 5, mu = 10, earlyStopping_max = 10, Sig_B = data$Sig_B, etol = 1e-4, noise = "additive", penalty = "lasso", proj = proj, constrain = constrain)
            }else {
                Zp <- ref_transform(data$Z)
                mu_u = rep(0,p)
                Sigma_x_lat = AR_covariance_matrix(p, rho)
                Sigma_w_lat = Sigma_x_lat+data$Sig_B
                supp_data = list(mu_x = theta, Sigma_x = Sigma_x_lat, Sigma_w = Sigma_w_lat, mu_w = theta)
                fit_additive <- fn_proposed(VV = Zp,y = data$y, alpha_real=beta_star, W = exp(data$Z),mu_u = mu_u, Sigma_u = data$Sig_B,EstimateSigma = T, Noestimate.data = supp_data, p_val = fpr)
                # print(fit_additive$beta.opt[1:10])
                # print(beta_star[1:10])
                if(model){
                    fit_additive$beta.opt = fit_additive$beta.lasso
                }
            }   
            # naive evaluation
            results_df$lambda[i] <- fit_additive$lambda.opt

            centered_X <- data$X - rep(fit_additive$mean.Z, each = n)
            # centered_X_ref <- ref_transform(centered_X)
            # measures <- eval_results(fit_additive$beta.opt[1:(p-1)], centered_X_ref, beta_star[1:(p-1)], beta_test = fit_additive$beta.test[1:(p-1)])
            measures <- eval_results(fit_additive$beta.opt, centered_X, beta_star, beta_test = fit_additive$beta.test)
            results_df$SE[i] <- measures$SE
            results_df$PE[i] <- measures$PE
            results_df$FPR[i] <- measures$FPR
            results_df$FNR[i] <- measures$FNR
            results_df$TPR[i] <- 1 - measures$FNR
            results_df$linf[i] <- measures$l_inf
            results_df$sum_beta[i] <- sum(fit_additive$beta.opt)

            results_bias[i,] <- fit_additive$beta.opt[c(1:10, p)]
        }
        runtime <- Sys.time() - start_time

        #--------------------------
        # bootstrap and p value

        bootstrap_mean <- apply(results_df, 2, mean)
        average_bias <- apply(results_bias, 2, mean)
        bootstrap_mean_std <- apply(results_df, 2, sd) / sqrt(N_sim)
        average_bias_std <- apply(results_bias, 2, sd) / sqrt(N_sim)
        p_value <- t.test(results_df$sum_beta, mu = 0)$p.value

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
        # write.table(sum_values,
        #     file = file_name_sum, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
        # )

        # print results
        evals <- rbind(bootstrap_mean, bootstrap_mean_std)
        colnames(evals) <- colnames(results_df)
        print(subdir_name)
        print(evals)

        print(runtime)
    }
}

