#include <Rcpp.h>
#include <RcppEigen.h>
#include <SymEigs.h>
#include <iostream>
#include <math.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/bessel.hpp>

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;
using namespace Spectra;
using namespace Eigen;
using namespace Rcpp;
using namespace std;
using boost::math::cyl_bessel_k;
typedef Map<MatrixXd> MapMatd;


void decomposeSymmetricMatrix(const Eigen::MatrixXd &M,
                              const int ncv,
                              const int k,
                              Eigen::VectorXd &lambda,
                              Eigen::MatrixXd &gamma) {
  /*
   * Parameters:
   *  - M: decomposed symmetric matrix
   *  - k:	Number of eigenvalues requested. 
   *     - Should satisfy 1≤nev≤n−1, where n is the size of matrix.
   *  - ncv: controls the convergence speed of the algorithm 
   * 
   */
  DenseSymMatProd<double> op(M);
  
  // Construct eigen solver object, requesting the largest k eigenvalues
  SymEigsSolver< double, LARGEST_ALGE, DenseSymMatProd<double>> eigs(&op, k, ncv);
  
  eigs.init();
  eigs.compute(1000, 1e-10);
  
  lambda = eigs.eigenvalues();
  gamma.noalias() = eigs.eigenvectors();
}

double thinPlateSplines(const double dist, const int d) {
  /*
   * Parameters:
   *  - dist: distance
   *  - d: dimension of positions (d <= 3)
   * Returns:
   *  - thin plate splines function of `dist`
   */
  double ret;
  
  if(d == 1)
    ret = pow(dist, 3) / 12;
  else if(d == 2) {
    if(dist != 0)
      ret = dist * dist * log(dist) / (8.0 * M_PI);
    else
      ret = 0;
  }
  else if(d == 3)
    ret = - dist / 8;
  else
    Rcpp::stop("Invalid dimension\n");
  
  return ret;
}

void createThinPlateMatrix(const MatrixXd s, MatrixXd &L) {
  /*
   * Parameters:
   *  - s: position matrix (n x d)
   *  - L: (return) resultant matrix (n x n)
   */
  int n(s.rows()), d(s.cols());
  double dist;
  // Update elements in the upper triangle 
  for(unsigned int i = 0; i < n; ++i) {
    for(unsigned int j = i + 1; j < n; ++j) {
      dist = (s.row(i) - s.row(j)).norm();
      L(i, j) = thinPlateSplines(dist, d);
    }
  }
  
  L += L.transpose().eval();
}

void predictThinPlateMatrix(const MatrixXd s_new,
                            const MatrixXd s,
                            MatrixXd &L) {
  /*
   * Parameters:
   *  - s_new: new position matrix (n1 x d)
   *  - s: reference position matrix (n2 x d)
   *  - L: (return) resultant matrix (n1 x n2)
   */
  int n1(s_new.rows()), n2(s.rows()), d(s.cols());
  double dist;
  
  for(unsigned int i = 0; i < n1; ++i) {
    for(unsigned int j = 0; j < n2; ++j) {
      dist = (s_new.row(i) - s.row(j)).norm();
      L(i, j) = thinPlateSplines(dist, d);
    }
  }
}

void mrtsBasis(const Eigen::MatrixXd s,
               const int k,
               Eigen::MatrixXd &Phi,
               Eigen::MatrixXd &B,
               Eigen::MatrixXd &BBB,
               Eigen::VectorXd &lambda,
               Eigen::MatrixXd &gamma) {
  /* MRTS: eigendecompose (I-B(B'B)^{-1}B')\Phi(I-B(B'B)^{-1}B')
   *    - Phi: thin plate splines
   *    - B := (1, s)
   *    
   * Parameters:
   *  - s: position matrix (n x d)
   *  - k:	Number of eigenvalues requested. 
   *     - Should satisfy 1≤nev≤n−1, where n is the size of matrix.
   * 
   */
  int n(s.rows()), d(s.cols());
  B = MatrixXd::Ones(n, d + 1);
  Phi = MatrixXd::Zero(n, n);
  
  // Create thin plate splines Phi
  createThinPlateMatrix(s, Phi);
  
  B.rightCols(d) = s;
  const Eigen::MatrixXd Bt = B.transpose();
  Eigen::MatrixXd BtB(MatrixXd(d + 1, d + 1).setZero().
                        selfadjointView<Lower>().rankUpdate(Bt));
  
  const Eigen::LLT<MatrixXd> llt(BtB);
  // BBB := B(B'B)^{-1}B'
  BBB = llt.solve(Bt);
  // Phi_proj := \Phi((I-B(B'B)^{-1}B')
  const Eigen::MatrixXd Phi_proj = Phi - (Phi * B) * BBB;
  // quadratic := ((I-B(B'B)^{-1}B')\Phi((I-B(B'B)^{-1}B')
  const Eigen::MatrixXd quadratic = Phi_proj - BBB.transpose() * (Bt * Phi_proj);

  // Set a convergence threshold for eigen-decomposition
  int ncv = min(n, max(2 * k + 1, 20));
  // Update lambda and gamma
  decomposeSymmetricMatrix(quadratic, ncv, k, lambda, gamma);
}

void mrtsCoreX(const Eigen::MatrixXd s,
              const Eigen::MatrixXd Phi,
              const Eigen::MatrixXd B,
              const Eigen::MatrixXd BBB,
              const Eigen::VectorXd lambda,
              const Eigen::MatrixXd gamma,
              int k,
              Eigen::MatrixXd &X,
              Eigen::VectorXd &nconst) {
  
  int n(s.rows()), d(s.cols());
  double root = sqrt(n);
  X = MatrixXd::Ones(n, k + d + 1);
  
  Eigen::MatrixXd X_center = s.rowwise() - s.colwise().mean();
  nconst = X_center.colwise().norm();
  
  X.block(0, 1, n, d) = X_center.array().rowwise() * (root / nconst.transpose().array());
  X.block(0, d + 1, n, k) = gamma * root;

  nconst /= root;
}

void mrtsCoreUZ(const Eigen::MatrixXd s,
              const Eigen::MatrixXd xobs_diag,
              const Eigen::MatrixXd Phi,
              const Eigen::MatrixXd B,
              const Eigen::MatrixXd BBB,
              const Eigen::VectorXd lambda,
              const Eigen::MatrixXd gamma,
              int k,
              Eigen::VectorXd &nconst) {
  
  int n(s.rows()), d(s.cols());
  double root = sqrt(n);
  Eigen::MatrixXd gammas;
  
  gammas = (gamma - B * (BBB * gamma)).array().rowwise() / lambda.transpose().array() * root;
  
  UZ = MatrixXd::Zero(n + d + 1, k + d + 1);
  UZ.block(0, 0, n, k) = gammas;
  UZ(n , k) = 1;
  UZ.block(n + 1, k + 1, d, d) = xobs_diag;
}



// [[Rcpp::export]]
Rcpp::List mrtsRcpp(const Eigen::Map<Eigen::MatrixXd> s,
                    const Eigen::Map<Eigen::MatrixXd> xobs_diag,
                    const int k) {
  Eigen::MatrixXd Phi, X, UZ, B, BBB, gamma;
  Eigen::VectorXd lambda, nconst;
  
  // Update B, BBB, lambda, gamma
  mrtsBasis(s, k, Phi, B, BBB, lambda, gamma);
  // Update X, UZ, nconst
  mrtsCoreX(s, Phi, B, BBB, lambda, gamma, k, X, nconst);
  mrtsCoreUZ(s, xobs_diag, Phi, B, BBB, lambda, gamma, k, UZ);
  
  return Rcpp::List::create(Rcpp::Named("X") = X,
                            Rcpp::Named("UZ") = UZ,
                            Rcpp::Named("BBBH") = BBB * Phi,
                            Rcpp::Named("nconst") = nconst);
}

// [[Rcpp::export]]
Rcpp::List predictMrtsRcpp(const Eigen::Map<Eigen::MatrixXd> s,
                           const Eigen::Map<Eigen::MatrixXd> xobs_diag,
                           const Eigen::Map<Eigen::MatrixXd> s_new,
                           const int k) {
  int n(s.rows()), d(s.cols()), n2(s_new.rows());
  Eigen::MatrixXd Phi, X, UZ, B, BBB, gamma, Phi_new;
  Eigen::VectorXd lambda, nconst;
  
  // Update B, BBB, lambda, gamma
  mrtsBasis(s, k, Phi, B, BBB, lambda, gamma);
  // Update X, UZ, nconst
  mrtsCoreX(s, Phi, B, BBB, lambda, gamma, k, X, nconst);
  mrtsCoreUZ(s, xobs_diag, Phi, B, BBB, lambda, gamma, k, UZ);
  
  Phi_new = MatrixXd::Zero(n2, n);
  //Create thin plate splines, Phi_new by new positions `s_new`
  predictThinPlateMatrix(s_new, s, Phi_new);
  
  Eigen::MatrixXd X1 = Phi_new * UZ.block(0, 0, n, k);
  Eigen::MatrixXd B_new = MatrixXd::Ones(n2, d + 1);
  B_new.rightCols(d) = s_new;
  
  return Rcpp::List::create(Rcpp::Named("X") = X,
                            Rcpp::Named("UZ") = UZ,
                            Rcpp::Named("BBBH") = BBB * Phi,
                            Rcpp::Named("nconst") = nconst,
                            Rcpp::Named("X1") = X1 - B_new * ((BBB * Phi) * UZ.block(0, 0, n, k)));
  
}

// [[Rcpp::export]]
Rcpp::List predictMrtsRcppWithBasis(const Eigen::Map<Eigen::MatrixXd> s,
                                    const Eigen::Map<Eigen::MatrixXd> xobs_diag,
                                    const Eigen::Map<Eigen::MatrixXd> s_new,
                                    const Eigen::Map<Eigen::MatrixXd> BBBH,
                                    const Eigen::Map<Eigen::MatrixXd> UZ,
                                    const Eigen::Map<Eigen::VectorXd> nconst,
                                    const int k) {
  int n(s.rows()), d(s.cols()), n2(s_new.rows());
  Eigen::MatrixXd Phi_new = MatrixXd::Zero(n2, n);
  //Create thin plate splines, Phi_new by new positions `s_new`
  predictThinPlateMatrix(s_new, s, Phi_new);
  
  Eigen::MatrixXd X1 = Phi_new * UZ.block(0, 0, n, k);
  Eigen::MatrixXd B = MatrixXd::Ones(n2, d + 1);
  B.rightCols(d) = s_new;
  
  return Rcpp::List::create(Rcpp::Named("X") = s,
                            Rcpp::Named("UZ") = UZ,
                            Rcpp::Named("BBBH") = BBBH,
                            Rcpp::Named("nconst") = nconst,
                            Rcpp::Named("X1") = X1 - B * ((BBBH) * UZ.block(0, 0, n, k)));
  
}

Eigen::MatrixXf maternCov(const Eigen::Map<Eigen::MatrixXd> s,
                          double tau,
                          double nu,
                          double rho) {
  /* Matern covariance function by Euclidean norm
   * 
   * Parameters
   *    - s: position matrix 
   *      - row: size
   *      - column: location dim)
   *    - tau: marginal variance (> 0)
   *    - nu: smoothness parameter (> 0)
   *    - rho: scale parameter (> 0)
   *
   * Returns
   *    covariance matrix (n x n)
   */
  int n(s.rows());
  double dist, scaled_dist, scalar, scaled_dist_nu, bessel;
  Eigen::MatrixXf cov(n, n);
  
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j <= i; ++j) {
      if (i == j)
        cov(i, j) = pow(tau, 2);
      else {
        dist = (s.row(i) - s.row(j)).norm();
        scaled_dist = (pow(2 * nu, 0.5) / rho) * dist;
        // Matern = scalar * (scaled_dist)^nu * bessel(scaled_dist)
        scalar = pow(tau, 2) * (pow(2, nu - 1) * std::tgamma(nu));
        scaled_dist_nu = pow(scaled_dist, nu);
        bessel = boost::math::cyl_bessel_k(nu, scaled_dist);
        cov(i, j) = scalar * scaled_dist_nu * bessel; 
      }
      cov(j, i) = cov(i, j);
    }
  }
  
  return cov;
}

// [[Rcpp::export]]
Eigen::MatrixXf inverseR(const Eigen::Map<Eigen::MatrixXd> s,
                         const double tau,
                         const double nu,
                         const double rho,
                         const double sigma2) {
  /* Inverse covariance matrix 
   *     - R = SNR * cov_matern + identity matrix
   *        - SNR: signal-to-noise ratio = rho^2 /sigma^2
   * Parameters
   *    - s: position matrix (n x d)
   *    - tau: marginal variance (> 0)
   *    - nu: smoothness parameter (> 0)
   *    - rho: scale parameter (> 0)
   *    - sigma2: noise variance (> 0)
   * Returns
   *    inverse covariance matrix (n x n)
   */
  int n(s.rows());
  Eigen::MatrixXf identity, R, R_inv, matern;
  identity = MatrixXf::Identity(n, n);
  
  matern = maternCov(s, tau, nu, rho);
  R = (pow(rho, 2) / sigma2) * matern + identity;
  R_inv = R.llt().solve(identity);
  
  return R_inv;
}


// [[Rcpp::export]]
double negLogLikelihood(const Eigen::Map<Eigen::MatrixXd> s,
                        double tau,
                        double nu,
                        double rho) {

  int n(s.rows()), d(s.cols());
  Eigen::MatrixXd Phi, F, UZ, B, BBB, gamma, Phi_new, R_inverse;
  Eigen::VectorXd lambda, nconst;
  
  // Update B, BBB, lambda, gamma
  mrtsBasis(s, k, Phi, B, BBB, lambda, gamma);
  // Fetch basis functions
  mrtsCoreX(s, Phi, B, BBB, lambda, gamma, k, F, nconst);
  // Fetch inverse of R
  R_inverse = inverseR(s, tau, nu, rho, sigma2);
  
}