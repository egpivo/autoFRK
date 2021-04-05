// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "autoFRK_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// eigenDecompose
Rcpp::List eigenDecompose(const Eigen::Map<Eigen::MatrixXd> matrix);
RcppExport SEXP _autoFRK_eigenDecompose(SEXP matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type matrix(matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenDecompose(matrix));
    return rcpp_result_gen;
END_RCPP
}
// mrtsrcpp
Rcpp::List mrtsrcpp(const Eigen::Map<Eigen::MatrixXd> Xu, const Eigen::Map<Eigen::MatrixXd> xobs_diag, const int k);
RcppExport SEXP _autoFRK_mrtsrcpp(SEXP XuSEXP, SEXP xobs_diagSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Xu(XuSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type xobs_diag(xobs_diagSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(mrtsrcpp(Xu, xobs_diag, k));
    return rcpp_result_gen;
END_RCPP
}
// mrtsrcpp_predict0
Rcpp::List mrtsrcpp_predict0(const Eigen::Map<Eigen::MatrixXd> Xu, const Eigen::Map<Eigen::MatrixXd> xobs_diag, const Eigen::Map<Eigen::MatrixXd> xnew, const int k);
RcppExport SEXP _autoFRK_mrtsrcpp_predict0(SEXP XuSEXP, SEXP xobs_diagSEXP, SEXP xnewSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Xu(XuSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type xobs_diag(xobs_diagSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type xnew(xnewSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(mrtsrcpp_predict0(Xu, xobs_diag, xnew, k));
    return rcpp_result_gen;
END_RCPP
}
// mrtsrcpp_predict
Rcpp::List mrtsrcpp_predict(const Eigen::Map<Eigen::MatrixXd> Xu, const Eigen::Map<Eigen::MatrixXd> xobs_diag, const Eigen::Map<Eigen::MatrixXd> xnew, const Eigen::Map<Eigen::MatrixXd> BBBH, const Eigen::Map<Eigen::MatrixXd> UZ, const Eigen::Map<Eigen::VectorXd> nconst, const int k);
RcppExport SEXP _autoFRK_mrtsrcpp_predict(SEXP XuSEXP, SEXP xobs_diagSEXP, SEXP xnewSEXP, SEXP BBBHSEXP, SEXP UZSEXP, SEXP nconstSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Xu(XuSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type xobs_diag(xobs_diagSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type xnew(xnewSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type BBBH(BBBHSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type UZ(UZSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type nconst(nconstSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(mrtsrcpp_predict(Xu, xobs_diag, xnew, BBBH, UZ, nconst, k));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_autoFRK_eigenDecompose", (DL_FUNC) &_autoFRK_eigenDecompose, 1},
    {"_autoFRK_mrtsrcpp", (DL_FUNC) &_autoFRK_mrtsrcpp, 3},
    {"_autoFRK_mrtsrcpp_predict0", (DL_FUNC) &_autoFRK_mrtsrcpp_predict0, 4},
    {"_autoFRK_mrtsrcpp_predict", (DL_FUNC) &_autoFRK_mrtsrcpp_predict, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_autoFRK(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
