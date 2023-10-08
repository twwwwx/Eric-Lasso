source("R/eic.R")
source("simulation/generate_data.R")
library("boot")
# library(devtools)
# install_github("zhaoni153/MicroBias")
library(limSolve)
library(MicroBias)
data(brooks)
names(brooks)
library(ggplot2)
library(reshape2)

set.seed(1234567)
model_list <- list()
model_list[["eico"]] <- c(TRUE, TRUE)
model_list[["coda"]] <- c(TRUE, FALSE)
model_list[["coco"]] <- c(FALSE, TRUE)
######### data preparation #########
delta <- brooks$delta
meta.data <- brooks$meta.data
obs.counts <- brooks$obs.counts
true.freq <- brooks$true.freq
Gene_names <- c("Avaginae", "Gvaginalis", "Lcrispatus", "Liners", "Pbivia", "Samnii", "GroupBStrep")

otu.tab <- obs.counts[, -which(colnames(obs.counts) == "Other")]
colnames(otu.tab) <- Gene_names
# Remove taxa that are not present by design
otu.tab <- otu.tab * delta

## Calculate relative abundancies
multi.otu <- which(rowSums(delta) > 1)
otu.tab <- otu.tab[multi.otu, ]
otu.tab <- otu.tab + 0.5
seq.length <- rowSums(otu.tab)
otu.tab <- otu.tab / rowSums(otu.tab)

## delete rows that only have one non-zero entry
true.freq <- true.freq[multi.otu, ]
meta.data <- meta.data[multi.otu, ]
delta <- delta[multi.otu, ]
k.s <- rowSums(delta)

## relative abundance to log matrix
Z <- log(as.matrix(otu.tab))
n <- nrow(Z)
p <- ncol(Z)
eps <- 0.5 / seq.length
X <- as.matrix(true.freq) + as.matrix(eps, n, 1) %*% matrix(1, 1, p)
X <- log(X)
clt <- function(X) {
    return(X - rowMeans(X))
}
Y <- clt(Z - X)

# x <- cbind(delta, model.matrix(~ factor(meta.data$Plate) + 0)[, 2:6])
# x <- true.freq
# mod <- MicroBias.fit(otu.tab = otu.tab, x = x, offset = true.freq, delta = delta)
# plot(seq(-1, -1, 1), seq(-1, -1, 1), xlab = "True", ylab = "Estimated", main = "EIC-Lasso", type = "l")


#################### EIC-Lasso ####################
N_bs <- 100
beta_matrix <- array(NA, dim = c(N_bs, 7, 7))
for (i in seq_along(model_list)) {
    model <- model_list[[i]]
    constrain <- model[1]
    proj <- model[2]
    model_name <- names(model_list)[i]
    begin.time <- Sys.time()
    for (n_bs in 1:N_bs) {
        if (n_bs %% 20 == 0) print(paste(n_bs, "%"))
        ## bootstrap
        bs.size <- round(n / 2)
        bs.size <- bs.size - bs.size %% 5 # make sure the size is a multiple of 5 for cross validation
        idx <- sample(1:n, bs.size, replace = TRUE)
        Z_bs <- Z[idx, ]
        nn <- nrow(Z_bs)
        X_bs <- X[idx, ]
        Sig_B <- var(Y[idx, ])

        for (i in 1:p) {
            y_bs <- as.matrix(Y[idx, i])
            fit <- eic(Z_bs, y_bs, nn, p, Sig_B = Sig_B, K = 5, noise = "additive", penalty = "lasso", proj = T, constrain = T, scale.Z = FALSE, scale.y = FALSE)
            beta_matrix[n_bs, , i] <- fit$beta
            centered_X <- X_bs - rep(fit$mean.Z, each = nn)
            estimated_y <- centered_X %*% fit$beta
        }
    }


    mean_beta <- apply(beta_matrix, c(2, 3), mean)
    beta_sgn <- sign(beta_matrix)
    beta_sgn <- apply(beta_sgn, c(2, 3), mean)
    beta_selected <- beta_matrix != 0
    beta_selected_freq <- apply(beta_selected, c(2, 3), sum) / N_bs
    end.time <- Sys.time()
    print(end.time - begin.time)


    ################ heatmap image ####################
    rownames(beta_selected_freq) <- rownames(Sig_B)
    colnames(beta_selected_freq) <- colnames(Sig_B)
    rownames(beta_sgn) <- rownames(Sig_B)
    colnames(beta_sgn) <- colnames(Sig_B)


    df <- melt(beta_selected_freq)
    df.sgn <- melt(beta_sgn)
    df$value[df$value < 0.7] <- NA
    df$sgn <- df.sgn$value > 0
    df$sgn[is.na(df$value)] <- NA
    colnames(df)[colnames(df) == "value"] <- "frequency"

    ggplot(df, aes(Var1, Var2)) +
        geom_tile(aes(fill = frequency), colour = "white") +
        # geom_text(aes(label = ifelse(value > 0.7, value, FALSE)), color = "black") +
        geom_text(aes(label = ifelse(sgn > 0, "+", "-")), color = "white", size = 8) +
        scale_fill_gradient(low = "#FF7256", high = "#8B3E2F", limits = c(0.7, 1), na.value = "white") +
        theme(axis.title = element_blank()) +
        coord_equal() +
        guides(fill = guide_colourbar(
            barwidth = 0.5,
            barheight = 20,
        ))
    ggsave(paste(model_name, "_heatmap_brooks.png"), plot = last_plot(), device = "png", dpi = 300)
    #################### fit image ####################
    Sig_B <- var(Y)
    n. <- n + 5 - n %% 5
    id <- sample(1:n, n. - n, replace = FALSE)

    Z_add <- rbind(Z, Z[id, ])
    X_add <- rbind(X, X[id, ])
    Y_add <- rbind(Y, Y[id, ])
    Predicted_y <- matrix(NA, n., p)
    for (i in 1:p) {
        y <- as.matrix(Y_add[, i])
        fit <- eic(Z_add, y, n., p, Sig_B = Sig_B, K = 5, noise = "additive", penalty = "lasso", proj = T, constrain = T, scale.Z = FALSE, scale.y = FALSE)
        intercept <- fit$mean.y - sum(fit$mean.Z * fit$beta)
        Predicted_y[, i] <- X_add %*% fit$beta + intercept
    }

    y.yp <- data.frame("true" = c(Y_add), "pred" = c(Predicted_y), group = factor(rep(Gene_names, each = n.)))
    # y.yp_filtered <- subset(y.yp, true != 0)

    ggplot(data = y.yp) +
        geom_point(aes(x = true, y = pred, color = group)) +
        geom_abline(intercept = 0, slope = 1, color = "#8B7D6B", linetype = "dashed", ) +
        labs(x = "True y", y = "Predicted y") +
        xlim(-2.85, 1.8) +
        ylim(-2.85, 1.8) +
        coord_equal()



    # Save the plot as a PNG file
    ggsave(paste(model_name, "fit_brooks.png"), plot = last_plot(), device = "png", dpi = 300)
}
