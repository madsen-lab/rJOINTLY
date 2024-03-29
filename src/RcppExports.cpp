// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// matDiMult
SEXP matDiMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, int n_cores);
RcppExport SEXP _JOINTLY_matDiMult(SEXP ASEXP, SEXP BSEXP, SEXP n_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(matDiMult(A, B, n_cores));
    return rcpp_result_gen;
END_RCPP
}
// matTriMult
SEXP matTriMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, Eigen::Map<Eigen::MatrixXd> C, int n_cores);
RcppExport SEXP _JOINTLY_matTriMult(SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP n_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type C(CSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(matTriMult(A, B, C, n_cores));
    return rcpp_result_gen;
END_RCPP
}
// matQuadMult
SEXP matQuadMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, Eigen::Map<Eigen::MatrixXd> C, Eigen::Map<Eigen::MatrixXd> D, int n_cores);
RcppExport SEXP _JOINTLY_matQuadMult(SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP DSEXP, SEXP n_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type C(CSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(matQuadMult(A, B, C, D, n_cores));
    return rcpp_result_gen;
END_RCPP
}
// cdist
NumericMatrix cdist(NumericMatrix x);
RcppExport SEXP _JOINTLY_cdist(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cdist(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_JOINTLY_matDiMult", (DL_FUNC) &_JOINTLY_matDiMult, 3},
    {"_JOINTLY_matTriMult", (DL_FUNC) &_JOINTLY_matTriMult, 4},
    {"_JOINTLY_matQuadMult", (DL_FUNC) &_JOINTLY_matQuadMult, 5},
    {"_JOINTLY_cdist", (DL_FUNC) &_JOINTLY_cdist, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_JOINTLY(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
