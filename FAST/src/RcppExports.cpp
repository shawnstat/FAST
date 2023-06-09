// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fast
List fast(NumericMatrix Xr, NumericMatrix Gr, int r, double lambda_1r, double lambda_2r, int n_iter, double converge);
RcppExport SEXP _FAST_fast(SEXP XrSEXP, SEXP GrSEXP, SEXP rSEXP, SEXP lambda_1rSEXP, SEXP lambda_2rSEXP, SEXP n_iterSEXP, SEXP convergeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Xr(XrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Gr(GrSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_1r(lambda_1rSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_2r(lambda_2rSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    Rcpp::traits::input_parameter< double >::type converge(convergeSEXP);
    rcpp_result_gen = Rcpp::wrap(fast(Xr, Gr, r, lambda_1r, lambda_2r, n_iter, converge));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FAST_fast", (DL_FUNC) &_FAST_fast, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_FAST(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
