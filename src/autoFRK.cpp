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
                              Eigen::VectorXd &rho,
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
  
  rho = eigs.eigenvalues();
  gamma.noalias() = eigs.eigenvectors();
}

void createThinPlateMatrix(const MatrixXd P, MatrixXd &L) {
  /*
   * Parameters:
   *  - P: position matrix (n x d)
   */
  int n(P.rows()), d(P.cols());
  double dist;
  // Update elements in the upper triangle 
  for(unsigned int i = 0; i < n; ++i) {
    for(unsigned int j = i + 1; j < n; ++j) {
      dist = (P.row(i) - P.row(j)).norm();
      if(d == 1)
        L(i, j) = pow(dist, 3) / 12;
      else if(d == 2) {
        if(dist != 0)
          L(i, j) = dist * dist * log(dist) / (8.0 * M_PI);
        else
          L(i, j) = 0;
      }
      else if(d == 3)
        L(i, j) = - dist / 8;
    }
  }
  
  L += L.transpose().eval();
}

void predictThinPlateMatrix(const MatrixXd P_new,
                            const MatrixXd P,
                            MatrixXd &L) {
  /*
   * Parameters:
   *  - P_new: predicted position matrix (n1 x d)
   *  - P: reference position matrix (n2 x d)
   */
  int n1(P_new.rows()), n2(P.rows()), d(P.cols());
  double dist;
  
  for(unsigned int i = 0; i < n1; ++i) {
    for(unsigned int j = 0; j < n2; ++j) {
      dist = (P_new.row(i) - P.row(j)).norm();
      if(d == 1)
        L(i, j) = pow(dist, 3) / 12;
      else if(d == 2) {
        if(dist != 0)
          L(i, j) = dist * dist * log(dist) / (8.0 * M_PI);
        else
          L(i, j) = 0;
      }
      else if(d == 3)
        L(i, j) = - dist / 8;
    }
  }
}

void mrtsBasis(const Eigen::MatrixXd Xu,
               const int k,
               Eigen::MatrixXd &Phi,
               Eigen::MatrixXd &B,
               Eigen::MatrixXd &BBB,
               Eigen::VectorXd &rho,
               Eigen::MatrixXd &gamma) {
  /* MRTS: eigendecompose (I-B(B'B)^{-1}B')\Phi(I-B(B'B)^{-1}B')
   *    - Phi: thin plate splines
   *    - B := (1, Xu)
   *    
   * Parameters:
   *  - Xu: position matrix (n x d)
   *  - k:	Number of eigenvalues requested. 
   *     - Should satisfy 1≤nev≤n−1, where n is the size of matrix.
   * 
   */
  int n(Xu.rows()), d(Xu.cols());
  B = MatrixXd::Ones(n, d + 1);
  Phi = MatrixXd::Zero(n, n);
  
  // Create thin plate splines Phi
  createThinPlateMatrix(Xu, Phi);
  
  B.rightCols(d) = Xu;
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
  
  // rho: eigenvalues; gamma: eigenvectors
  int ncv = min(n, max(2 * k + 1, 20));
  decomposeSymmetricMatrix(quadratic, ncv, k, rho, gamma);
}

void mrtsCore(const Eigen::MatrixXd Xu,
              const Eigen::MatrixXd xobs_diag,
              const Eigen::MatrixXd Phi,
              const Eigen::MatrixXd B,
              const Eigen::MatrixXd BBB,
              const Eigen::VectorXd rho,
              const Eigen::MatrixXd gamma,
              int k,
              Eigen::MatrixXd &X,
              Eigen::MatrixXd &UZ,
              Eigen::VectorXd &nconst) {
  
  int n(Xu.rows()), d(Xu.cols());
  double root = sqrt(n);
  Eigen::MatrixXd gammas;
  X = MatrixXd::Ones(n, k + d + 1);
  
  gammas = (gamma - B * (BBB * gamma)).array().rowwise() / rho.transpose().array() * root;
  Eigen::MatrixXd X_center = Xu.rowwise() - Xu.colwise().mean();
  nconst = X_center.colwise().norm();
  
  X.block(0, 1, n, d) = X_center.array().rowwise() * (root / nconst.transpose().array());
  X.block(0, d + 1, n, k) = gamma * root;
  
  UZ = MatrixXd::Zero(n + d + 1, k + d + 1);
  UZ.block(0, 0, n, k) = gammas;
  UZ(n , k) = 1;
  UZ.block(n + 1, k + 1, d, d) = xobs_diag;
  nconst /= root;
}

// [[Rcpp::export]]
Rcpp::List mrtsRcpp(const Eigen::Map<Eigen::MatrixXd> Xu,
                    const Eigen::Map<Eigen::MatrixXd> xobs_diag,
                    const int k) {
  Eigen::MatrixXd Phi, X, UZ, B, BBB, gamma;
  Eigen::VectorXd rho, nconst;
  
  mrtsBasis(Xu, k, Phi, B, BBB, rho, gamma);
  mrtsCore(Xu, xobs_diag, Phi, B, BBB, rho, gamma, k, X, UZ, nconst);
  
  return Rcpp::List::create(Rcpp::Named("X") = X,
                            Rcpp::Named("UZ") = UZ,
                            Rcpp::Named("BBBH") = BBB * Phi,
                            Rcpp::Named("nconst") = nconst);
}

// [[Rcpp::export]]
Rcpp::List predictMrtsRcpp(const Eigen::Map<Eigen::MatrixXd> Xu,
                           const Eigen::Map<Eigen::MatrixXd> xobs_diag,
                           const Eigen::Map<Eigen::MatrixXd> x_new,
                           const int k) {
  int n(Xu.rows()), d(Xu.cols()), n2(x_new.rows());
  Eigen::MatrixXd Phi, X, UZ, B, BBB, gamma, Phi_new;
  Eigen::VectorXd rho, nconst;
  
  mrtsBasis(Xu, k, Phi, B, BBB, rho, gamma);
  mrtsCore(Xu, xobs_diag, Phi, B, BBB, rho, gamma, k, X, UZ, nconst);
  
  Phi_new = MatrixXd::Zero(n2, n);
  predictThinPlateMatrix(x_new, Xu, Phi_new);
  
  Eigen::MatrixXd X1 = Phi_new * UZ.block(0, 0, n, k);
  Eigen::MatrixXd B_new = MatrixXd::Ones(n2, d + 1);
  B_new.rightCols(d) = x_new;
  
  return Rcpp::List::create(Rcpp::Named("X") = X,
                            Rcpp::Named("UZ") = UZ,
                            Rcpp::Named("BBBH") = BBB * Phi,
                            Rcpp::Named("nconst") = nconst,
                            Rcpp::Named("X1") = X1 - B_new * ((BBB * Phi) * UZ.block(0, 0, n, k)));
  
}

// [[Rcpp::export]]
Rcpp::List predictMrtsRcppWithBasis(const Eigen::Map<Eigen::MatrixXd> Xu,
                                    const Eigen::Map<Eigen::MatrixXd> xobs_diag,
                                    const Eigen::Map<Eigen::MatrixXd> x_new,
                                    const Eigen::Map<Eigen::MatrixXd> BBBH,
                                    const Eigen::Map<Eigen::MatrixXd> UZ,
                                    const Eigen::Map<Eigen::VectorXd> nconst,
                                    const int k) {
  int n(Xu.rows()), d(Xu.cols()), n2(x_new.rows());
  Eigen::MatrixXd Phi_new = MatrixXd::Zero(n2, n);
  //Predict Tinn Plate Matrix by `x_new`
  predictThinPlateMatrix(x_new, Xu, Phi_new);
  
  Eigen::MatrixXd X1 = Phi_new * UZ.block(0, 0, n, k);
  Eigen::MatrixXd B = MatrixXd::Ones(n2, d + 1);
  B.rightCols(d) = x_new;
  
  return Rcpp::List::create(Rcpp::Named("X") = Xu,
                            Rcpp::Named("UZ") = UZ,
                            Rcpp::Named("BBBH") = BBBH,
                            Rcpp::Named("nconst") = nconst,
                            Rcpp::Named("X1") = X1 - B * ((BBBH) * UZ.block(0, 0, n, k)));
  
}

// [[Rcpp::export]]
Eigen::MatrixXf maternCov(const Eigen::Map<Eigen::MatrixXd> s,
                          double tau,
                          double nu,
                          double rho) {
  /* Matern covariance function by Euclidean norm
   * 
   * Parameters
   *    s: location matrix 
   *      - row: size
   *      - column: location dim)
   *    tau, nu, rho: positive real number
   * Returns
   *    covariance matrix (n x n)
   */
  int n(s.rows());
  double l2_dist, scalar, first, second, third;
  Eigen::MatrixXf cov(n, n);
  
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j <= i; ++j) {
      if (i == j) {
        cov(i, j) = pow(tau, 2);
      }
      else {
        l2_dist = (s.row(i) - s.row(j)).norm();
        scalar = (pow(2 * nu, 0.5) / rho) * l2_dist;
        first = pow(tau, 2) * (pow(2, nu - 1) * std::tgamma(nu));
        second = pow(scalar, nu);
        third = boost::math::cyl_bessel_k(nu, scalar);
        cov(i, j) = first * second * third; 
      }
      cov(j, i) = cov(i, j);
    }
  }
  
  return cov;
}
