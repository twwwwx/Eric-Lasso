#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Soft threshold function
double soft_threshold(double num, double lambda) {
  return std::copysign(std::max(std::abs(num) - lambda, 0.0), num);
}

// [[Rcpp::export]]
List lasso_covariance_con_cpp(int n,
                          int p,
                          double lambda,
                          arma::mat XX,
                          arma::vec Xy,
                          arma::vec beta_start,
                          double mu = 1,
                          std::string penalty = "lasso") {
  arma::vec beta = beta_start;
  arma::mat wp = beta;
  int m = 1;
  double alpha = 0;
  std::vector<double> loss_list;

  // Control parameters
  // int maxIter = control["maxIter"];
  // double optTol = control["optTol"];
  // double zeroThreshold = control["zeroThreshold"];
  int maxIter = 1000;
  double optTol = 1e-5;
  double zeroThreshold = 1e-6;

  // Compute the product of XX with beta
  arma::vec s = XX * beta;
  while (m < maxIter) {
    arma::vec beta_old = beta;
    for (int j = 0; j < p; j++) {
      // Compute the Shoot and Update the variable
      double S0 = s[j] - XX(j, j) * beta_old[j] - Xy[j];
      S0 = S0 + mu * (arma::sum(beta) - beta[j] + alpha);
      if (std::isnan(S0)) {
        beta[j] = 0;
        continue;
      }

      beta[j] = soft_threshold(-S0, lambda) / (XX(j, j) + mu);
      s = s + XX.col(j) * (beta[j] - beta_old[j]); // s: XX %*% beta_new
    }
    // Update
    wp.insert_cols(wp.n_cols, beta);
    alpha += arma::sum(beta);
    // Check termination for early stopping
    if (arma::sum(arma::abs(beta - beta_old)) < optTol) {
      break;
    }
    // Check the loss function
    double loss = 0.5 * arma::as_scalar(beta.t() * XX * beta) - arma::as_scalar(Xy.t() * beta) + lambda * arma::sum(arma::abs(beta)) + mu * std::pow(arma::sum(beta), 2);
    loss_list.push_back(loss);

    m++;
  }
  arma::vec w = beta;
  // We impose very small coefficients to be equal to zero
  w.elem(arma::find(arma::abs(w) < zeroThreshold)).zeros();
  return List::create(Named("coefficients") = w, Named("coef_list") = wp, Named("num_it") = m, Named("loss") = loss_list);
}