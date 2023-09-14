soft_threshold <- function(num, lambda) {
  sign_num <- sign(num)
  sign_num * max((abs(num) - lambda), 0)
}



#' Compositional Lasso in covariance form
#'
#' Solve the least squares loss with lasso penalty written in a form with the covariance matrix : \eqn{\frac{1}{2} \beta^{'} \Sigma \beta - \rho^{'} \beta + \lambda \|\beta\|_1} subject to sum(beta_j) = 0
#'
#' @param n Number of samples of the design matrix
#' @param p Number of features of the matrix
#' @param lambda penalty parameter
#' @param control Including control parameters : max of iterations, tolerance for the convergence of the error, zero threshold to put to zero small beta coefficients
#' @param XX Design matrix corresponding to \eqn{\frac{1}{n} X'X} or a modified version in the case of CoCoLasso
#' @param Xy Rho parameter corresponding to \eqn{\frac{1}{n} X'y} or a modified version in the case of CoCoLasso
#' @param beta.start Initial value of beta
#' @param penalty Type of penalty used : can be lasso penalty or SCAD penalty
#'
#' @return list containing \itemize{
#' \item coefficients : Coefficients corresponding to final beta after convergence of the algoritm
#' \item coef.list : Matrix of coefficients for beta for all iterations
#' \item num.it Number of iterations
#' }
#'
#' @export

lasso_covariance_con <- function(n,
                                 p,
                                 lambda,
                                 control = list(
                                   maxIter = 1000,
                                   optTol = 10^(-5),
                                   zeroThreshold = 10^(-6)
                                 ),
                                 XX,
                                 Xy,
                                 beta.start,
                                 mu = 1,
                                 penalty = c("lasso", "SCAD")) {
  beta <- beta.start
  wp <- beta
  m <- 1
  alpha <- 0
  loss_list <- list()

  ## compute the product of XX with beta
  s <- XX %*% beta
  while (m < control$maxIter) {
    beta_old <- beta
    for (j in 1:p) {
      # Compute the Shoot and Update the variable
      S0 <- s[j] - XX[j, j] * beta_old[j] - Xy[j]
      S0 <- S0 + mu * (sum(beta) - beta[j] + alpha)
      if (sum(is.na(S0)) >= 1) {
        beta[j] <- 0
        next
      }

      beta[j] <- soft_threshold(-S0, lambda) / (XX[j, j] + mu)
      s <- s + XX[, j] * (beta[j] - beta_old[j]) # s: XX %*% beta_new
    }
    # Update
    wp <- cbind(wp, beta)
    alpha <- alpha + sum(beta)
    # Check termination for early stopping
    if (sum(abs(beta - beta_old), na.rm = TRUE) < control$optTol) {
      break
    }
    # Check the loss function
    loss <- 0.5 * t(beta) %*% XX %*% beta - t(Xy) %*% beta + lambda * sum(abs(beta)) + mu * sum(beta)^2
    loss_list <- c(loss_list, loss)

    m <- m + 1
  }
  w <- beta
  # We impose very small coefficients to be equal to zero
  w[abs(w) < control$zeroThreshold] <- 0
  return(list(coefficients = w, coef.list = wp, num.it = m, loss = loss_list))
}
