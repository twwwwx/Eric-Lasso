#' Projected covariance matrices
#'
#' Creating for the cross validation step of the CoCoLasso
#'
#' @param K Number of folds for the cross validation
#' @param mat Covariance matrix to be projected
#' @param y Response vector
#' @param p Number of predictors
#' @param mu Penalty parameter for the ADMM algorithm.
#' @param ratio_matrix Observation matrix used in the missing data setting
#' @param etol Tolerance used in the ADMM algorithm
#' @param noise Type of setting chosen : additive or missing
#' @param mode ADMM or HM

#'
#' @return list containing \itemize{
#' \item \code{sigma_global} projected matrix for \code{mat}
#' \item \code{rho_global} rho parameter for \code{mat}
#' \item \code{list_matrices_lasso} list of the projected matrices for \code{mat} deprived of the k-fold during cross validation
#' \item \code{list_matrices_error} list of the projected matrices for the k-fold of \code{mat}
#' \item \code{list_rho_lasso} list of the modified \code{rho} for \code{mat} deprived of the K-fold during cross validation
#' \item \code{list_rho_error} list of the modified \code{rho} for the k-fold of \code{mat}
#' }
#'
#'
#' @export

cv_matrices_noproj <- function(K,
                               mat,
                               y,
                               p,
                               mu,
                               ratio_matrix = NULL,
                               etol = 1e-4,
                               noise = c("additive", "missing"),
                               mode = "ADMM") {
  # calculate the K nearest PSD covariance matrices
  #      in the cross validation process

  # there is no correlation between the error for different features.
  # @param ratio_matrix is the observation matrix. Its j,k term calculates the number of rows where we observe j and
  # k feature at the same time if j!=k, and the number of rows where we observe the j feature if j==k

  n <- nrow(mat)
  p <- ncol(mat)
  n_without_fold <- n - floor(n / K)
  n_one_fold <- floor(n / K)

  folds <- sample(cut(seq(1, n), breaks = K, labels = FALSE))
  list_matrices_lasso <- list()
  list_matrices_error <- list()
  list_rho_lasso <- list()
  list_rho_error <- list()

  ### Case where there is additive noise

  # We calculate the global nearest PSD cov matrix and the global surrogate rho, when we take into account the whole data set

  sigma_global <- 1 / n * t(mat) %*% mat
  rho_global <- t(mat) %*% y / n

  for (i in 1:K) {
    # We calculate the necessary matrices for the cross validation

    # print(paste("Doing the", i, "fold"))
    index <- which(folds == i, arr.ind = TRUE)

    mat_train <- mat[-index, ]
    mat_cov_train <- 1 / n_without_fold * t(mat_train) %*% mat_train

    list_matrices_lasso <- rlist::list.append(list_matrices_lasso, mat_cov_train)

    mat_test <- mat[index, ]
    mat_cov_test <- 1 / n_one_fold * t(mat_test) %*% mat_test

    list_matrices_error <- rlist::list.append(list_matrices_error, mat_cov_test)

    # Calculating the surrogate rho when we remove the kth fold, to resolve lasso problem during cross validation
    y_train <- y[-index, ]
    rho_train <- 1 / n_without_fold * t(mat_train) %*% y_train
    list_rho_lasso <- rlist::list.append(list_rho_lasso, rho_train)

    # Calculating the surrogate rho for the kth fold, to calculate the error on the problem solved without the kth fold
    y_test <- y[index, ]
    rho_test <- 1 / n_one_fold * t(mat_test) %*% y_test
    list_rho_error <- rlist::list.append(list_rho_error, rho_test)
  }


  return(list(
    sigma_global = sigma_global, rho_global = rho_global, list_matrices_lasso = list_matrices_lasso, list_matrices_error = list_matrices_error,
    list_rho_lasso = list_rho_lasso, list_rho_error = list_rho_error
  ))
}
