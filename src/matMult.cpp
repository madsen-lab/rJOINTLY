// [[Rcpp::depends(Rcpp, RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppEigen.h>
#include <Rcpp.h>
#include <string>

using namespace Rcpp;

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

// [[Rcpp::export]]
NumericMatrix cdist(NumericMatrix x){
  int n=x.nrow(),ncol=x.ncol(),i,j,k;
  NumericMatrix out(n,n);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++){
      double sum=0;
      for(k=0;k<ncol;k++)sum+=pow(x(i,k)-x(j,k),2);
      out(i,j)=sqrt(sum);
    }
    out;
}