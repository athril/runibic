// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// runibic_params
void runibic_params(double t, double q, double f, int nbic, int div);
RcppExport SEXP _runibic_runibic_params(SEXP tSEXP, SEXP qSEXP, SEXP fSEXP, SEXP nbicSEXP, SEXP divSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type f(fSEXP);
    Rcpp::traits::input_parameter< int >::type nbic(nbicSEXP);
    Rcpp::traits::input_parameter< int >::type div(divSEXP);
    runibic_params(t, q, f, nbic, div);
    return R_NilValue;
END_RCPP
}
// discretize
Rcpp::IntegerMatrix discretize(Rcpp::NumericMatrix x);
RcppExport SEXP _runibic_discretize(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(discretize(x));
    return rcpp_result_gen;
END_RCPP
}
// unisort
Rcpp::IntegerMatrix unisort(Rcpp::IntegerMatrix x);
RcppExport SEXP _runibic_unisort(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(unisort(x));
    return rcpp_result_gen;
END_RCPP
}
// pairwiseLCS
Rcpp::IntegerMatrix pairwiseLCS(Rcpp::IntegerVector x, Rcpp::IntegerVector y);
RcppExport SEXP _runibic_pairwiseLCS(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(pairwiseLCS(x, y));
    return rcpp_result_gen;
END_RCPP
}
// backtrackLCS
Rcpp::IntegerVector backtrackLCS(Rcpp::IntegerVector x, Rcpp::IntegerVector y);
RcppExport SEXP _runibic_backtrackLCS(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(backtrackLCS(x, y));
    return rcpp_result_gen;
END_RCPP
}
// calculateLCS
Rcpp::List calculateLCS(Rcpp::IntegerMatrix discreteInput, bool useFibHeap);
RcppExport SEXP _runibic_calculateLCS(SEXP discreteInputSEXP, SEXP useFibHeapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type discreteInput(discreteInputSEXP);
    Rcpp::traits::input_parameter< bool >::type useFibHeap(useFibHeapSEXP);
    rcpp_result_gen = Rcpp::wrap(calculateLCS(discreteInput, useFibHeap));
    return rcpp_result_gen;
END_RCPP
}
// cluster
Rcpp::List cluster(Rcpp::IntegerMatrix discreteInput, Rcpp::IntegerMatrix discreteInputValues, Rcpp::IntegerVector scores, Rcpp::IntegerVector geneOne, Rcpp::IntegerVector geneTwo, int rowNumber, int colNumber);
RcppExport SEXP _runibic_cluster(SEXP discreteInputSEXP, SEXP discreteInputValuesSEXP, SEXP scoresSEXP, SEXP geneOneSEXP, SEXP geneTwoSEXP, SEXP rowNumberSEXP, SEXP colNumberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type discreteInput(discreteInputSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type discreteInputValues(discreteInputValuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type scores(scoresSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type geneOne(geneOneSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type geneTwo(geneTwoSEXP);
    Rcpp::traits::input_parameter< int >::type rowNumber(rowNumberSEXP);
    Rcpp::traits::input_parameter< int >::type colNumber(colNumberSEXP);
    rcpp_result_gen = Rcpp::wrap(cluster(discreteInput, discreteInputValues, scores, geneOne, geneTwo, rowNumber, colNumber));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_runibic_runibic_params", (DL_FUNC) &_runibic_runibic_params, 5},
    {"_runibic_discretize", (DL_FUNC) &_runibic_discretize, 1},
    {"_runibic_unisort", (DL_FUNC) &_runibic_unisort, 1},
    {"_runibic_pairwiseLCS", (DL_FUNC) &_runibic_pairwiseLCS, 2},
    {"_runibic_backtrackLCS", (DL_FUNC) &_runibic_backtrackLCS, 2},
    {"_runibic_calculateLCS", (DL_FUNC) &_runibic_calculateLCS, 2},
    {"_runibic_cluster", (DL_FUNC) &_runibic_cluster, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_runibic(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
