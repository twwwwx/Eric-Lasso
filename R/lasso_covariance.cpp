#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List lasso_covariance_rcpp(int n, int p, double lambda, 
                            arma::mat XX, arma::vec Xy, arma::vec beta_start, 
                            std::string penalty = "lasso") {
    arma::vec beta = beta_start;
    arma::mat wp = beta;
    int m = 1;
    
    arma::vec s = XX * beta;
    double lambda0 = lambda;
    
    // int maxIter = control["maxIter"],
    //     optTol = control["optTol"],
    //     zeroThreshold = control["zeroThreshold"];
    int maxIter = 1000;
    double optTol = 1e-5;
    double zeroThreshold = 1e-6;
    while (m < maxIter) {
        arma::vec beta_old = beta;
        
        for (int j = 0; j < p; j++) {
            double S0 = s[j] - XX(j, j) * beta_old[j] - Xy[j];
            if (std::isnan(S0)) {
                beta[j] = 0;
                continue;
            }
            
            double w_j = 1.0;
            if (penalty == "SCAD") {
                double a = 3.7;
                if (std::abs(beta[j]) > lambda && std::abs(beta[j]) <= a * lambda) {
                    w_j = (a * lambda - std::abs(beta[j])) / (lambda * (a - 1));
                } else if (std::abs(beta[j]) > a * lambda) {
                    w_j = 0.0;
                }
            }
            
            lambda = w_j * lambda0;
            if (S0 > lambda) {
                beta[j] = (lambda - S0) / XX(j, j);
                s += XX.col(j) * (beta[j] - beta_old[j]);
            } else if (S0 < -lambda) {
                beta[j] = (-lambda - S0) / XX(j, j);
                s += XX.col(j) * (beta[j] - beta_old[j]);
            } else {
                beta[j] = 0;
            }
        }
        
        wp.insert_cols(wp.n_cols, beta);
        if (arma::sum(arma::abs(beta - beta_old)) < optTol) {
            break;
        }
        m++;
    }
    
    beta.elem(find(arma::abs(beta) < zeroThreshold)).zeros();
    
    return List::create(
        Named("coefficients") = beta,
        Named("coef.list") = wp,
        Named("num.it") = m
    );
}
