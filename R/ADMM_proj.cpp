#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Efficient projection onto L1 ball
// [[Rcpp::export]]
arma::vec l1proj_cpp(const arma::vec& v, double b) {
    // Ensure radius b is positive
    if (b <= 0) {
      stop("b must be greater than 0");
    }
  
    // Sort the absolute values of v in decreasing order
    arma::vec u = sort(abs(v), "descend");
  
    // Compute the cumulative sum of sorted values
    arma::vec sv = cumsum(u);
  
    // Find the maximum rho such that u[rho] > (sv[rho] - b) / (rho + 1)
    arma::uword rho = 0;
    for (arma::uword i = 0; i < u.n_elem; ++i) {
      if (u[i] > (sv[i] - b) / (i + 1)) {
        rho = i;
      } else {
        break;
      }
    }
  
    // Compute the threshold theta
    double theta = std::max(0.0, (sv[rho] - b) / (rho + 1));
  
    // Compute the projection w
    arma::vec w = arma::sign(v) % arma::max(arma::abs(v) - theta, arma::zeros<arma::vec>(v.n_elem));
  
    return w;
  }

// ADMM algorithm for nearest positive semi-definite matrix
// [[Rcpp::export]]
List ADMM_proj_cpp(arma::mat mat, double epsilon = 1e-4, double mu = 10, int it_max = 1000,
               double etol = 1e-4, double etol_distance = 1e-4) {
    int p = mat.n_rows;
    arma::mat R = diagmat(mat);
    arma::mat S = arma::zeros(p, p);
    arma::mat L = arma::zeros(p, p);
    
    int itr = 0;
    arma::vec iteration, eps_R, eps_S, eps_primal, time, distance;
    
    while (itr < it_max) {
        arma::mat Rp = R;
        arma::mat Sp = S;
        // auto start = std::chrono::steady_clock::now();
        
        // Print matrix dimensions
        // std::cout << "Iteration " << itr << ":\n";
        // std::cout << "  S: " << arma::size(S) << "\n";
        // std::cout << "  L: " << arma::size(L) << "\n";
        
        // R step
        arma::mat W = mat + S + mu * L;
        arma::vec eigval;
        arma::mat eigvec;
        eig_sym(eigval, eigvec, W);
        // std::cout << "  eigval: " << arma::size(eigval) << "\n";
        // std::cout << "  eigvec: " << arma::size(eigvec) << "\n";
        
        // arma::mat deigval = diagmat(eigval.transform([epsilon](double x) { return std::max(x, epsilon); }));
        // std::cout << "  deigval: " << arma::size(deigval) << "\n";
        // std::cout << "  R: " << arma::size(R) << "\n";
        R = eigvec * diagmat(eigval.transform([epsilon](double x) { return std::max(x, epsilon); })) * eigvec.t();
        
        
        // S step
        arma::mat M = R - mat - mu * L;
        // arma::vec M_lower = M.elem(find(trimatu(ones(p, p), 0)));
        arma::uvec lower_indices = arma::trimatl_ind(size(M), -1);
        arma::vec M_lower = M(lower_indices);
        // std::cout << "  M_lower: " << arma::size(M_lower) << "\n";
        arma::vec S_lower = M_lower - l1proj_cpp(M_lower, mu / 2);
        // std::cout << "  S_lower: " << arma::size(S_lower) << "\n";
        S(lower_indices) = S_lower;
        S = symmatu(S);
        
        // L step
        // std::cout << "  S: " << arma::size(S) << "\n";
        // std::cout << "  tmp: " << arma::size(tmp) << "\n";
        L -= (R - S - mat) / mu;

        if ((abs(R - Rp).max() < etol && abs(S - Sp).max() < etol && abs(R - S - mat).max() < etol) || 
        (std::abs((Rp - mat).max() - (R - mat).max()) < etol_distance)) {
            break;
        }
        
    
        itr++;
        if (itr % 20 == 0) {
            mu /= 2;
        }
    }
    
    return List::create(Named("mat") = R);
}