// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP matDiMult(const Eigen::Map<Eigen::MatrixXd> A,
               Eigen::Map<Eigen::MatrixXd> B,
               int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP matTriMult(const Eigen::Map<Eigen::MatrixXd> A,
                Eigen::Map<Eigen::MatrixXd> B,
                Eigen::Map<Eigen::MatrixXd> C,
                int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd D = A * B * C;
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP matQuadMult(const Eigen::Map<Eigen::MatrixXd> A,
                Eigen::Map<Eigen::MatrixXd> B,
                Eigen::Map<Eigen::MatrixXd> C,
                Eigen::Map<Eigen::MatrixXd> D,
                int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd E = A * B * C * D;
  return Rcpp::wrap(E);
}