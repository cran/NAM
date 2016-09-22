// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// calcSize
int calcSize(NumericVector col, NumericVector fam);
RcppExport SEXP NAM_calcSize(SEXP colSEXP, SEXP famSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type col(colSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fam(famSEXP);
    rcpp_result_gen = Rcpp::wrap(calcSize(col, fam));
    return rcpp_result_gen;
END_RCPP
}
// emBA
SEXP emBA(NumericVector y, NumericMatrix gen, double df, double R2, int it);
RcppExport SEXP NAM_emBA(SEXP ySEXP, SEXP genSEXP, SEXP dfSEXP, SEXP R2SEXP, SEXP itSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    rcpp_result_gen = Rcpp::wrap(emBA(y, gen, df, R2, it));
    return rcpp_result_gen;
END_RCPP
}
// emBB
SEXP emBB(NumericVector y, NumericMatrix gen, double df, double R2, int it, double Pi);
RcppExport SEXP NAM_emBB(SEXP ySEXP, SEXP genSEXP, SEXP dfSEXP, SEXP R2SEXP, SEXP itSEXP, SEXP PiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type Pi(PiSEXP);
    rcpp_result_gen = Rcpp::wrap(emBB(y, gen, df, R2, it, Pi));
    return rcpp_result_gen;
END_RCPP
}
// emBC
SEXP emBC(NumericVector y, NumericMatrix gen, double df, double R2, int it, double Pi);
RcppExport SEXP NAM_emBC(SEXP ySEXP, SEXP genSEXP, SEXP dfSEXP, SEXP R2SEXP, SEXP itSEXP, SEXP PiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type Pi(PiSEXP);
    rcpp_result_gen = Rcpp::wrap(emBC(y, gen, df, R2, it, Pi));
    return rcpp_result_gen;
END_RCPP
}
// emRR
SEXP emRR(NumericVector y, NumericMatrix gen, double df, double R2, int it);
RcppExport SEXP NAM_emRR(SEXP ySEXP, SEXP genSEXP, SEXP dfSEXP, SEXP R2SEXP, SEXP itSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    rcpp_result_gen = Rcpp::wrap(emRR(y, gen, df, R2, it));
    return rcpp_result_gen;
END_RCPP
}
// funI
NumericVector funI(NumericVector col, int fam, int finsiz, int f);
RcppExport SEXP NAM_funI(SEXP colSEXP, SEXP famSEXP, SEXP finsizSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type col(colSEXP);
    Rcpp::traits::input_parameter< int >::type fam(famSEXP);
    Rcpp::traits::input_parameter< int >::type finsiz(finsizSEXP);
    Rcpp::traits::input_parameter< int >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(funI(col, fam, finsiz, f));
    return rcpp_result_gen;
END_RCPP
}
// funX
NumericVector funX(NumericVector col, int finsiz);
RcppExport SEXP NAM_funX(SEXP colSEXP, SEXP finsizSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type col(colSEXP);
    Rcpp::traits::input_parameter< int >::type finsiz(finsizSEXP);
    rcpp_result_gen = Rcpp::wrap(funX(col, finsiz));
    return rcpp_result_gen;
END_RCPP
}
// gs
void gs(NumericMatrix C, NumericVector g, NumericVector r, int N);
RcppExport SEXP NAM_gs(SEXP CSEXP, SEXP gSEXP, SEXP rSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type C(CSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g(gSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    gs(C, g, r, N);
    return R_NilValue;
END_RCPP
}
// inputRow
NumericVector inputRow(NumericVector x, int exp, int n);
RcppExport SEXP NAM_inputRow(SEXP xSEXP, SEXP expSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type exp(expSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(inputRow(x, exp, n));
    return rcpp_result_gen;
END_RCPP
}
// KMUP
SEXP KMUP(NumericMatrix X, NumericVector b, NumericVector d, NumericVector xx, NumericVector E, NumericVector L, double Ve, double pi);
RcppExport SEXP NAM_KMUP(SEXP XSEXP, SEXP bSEXP, SEXP dSEXP, SEXP xxSEXP, SEXP ESEXP, SEXP LSEXP, SEXP VeSEXP, SEXP piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type Ve(VeSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    rcpp_result_gen = Rcpp::wrap(KMUP(X, b, d, xx, E, L, Ve, pi));
    return rcpp_result_gen;
END_RCPP
}
// KMUP2
SEXP KMUP2(NumericMatrix X, NumericVector Use, NumericVector b, NumericVector d, NumericVector xx, NumericVector E, NumericVector L, double Ve, double pi);
RcppExport SEXP NAM_KMUP2(SEXP XSEXP, SEXP UseSEXP, SEXP bSEXP, SEXP dSEXP, SEXP xxSEXP, SEXP ESEXP, SEXP LSEXP, SEXP VeSEXP, SEXP piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Use(UseSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type Ve(VeSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    rcpp_result_gen = Rcpp::wrap(KMUP2(X, Use, b, d, xx, E, L, Ve, pi));
    return rcpp_result_gen;
END_RCPP
}
// SAMP
void SAMP(NumericMatrix C, NumericVector g, NumericVector r, int N, double Ve);
RcppExport SEXP NAM_SAMP(SEXP CSEXP, SEXP gSEXP, SEXP rSEXP, SEXP NSEXP, SEXP VeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type C(CSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g(gSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type Ve(VeSEXP);
    SAMP(C, g, r, N, Ve);
    return R_NilValue;
END_RCPP
}
// SAMP2
void SAMP2(NumericMatrix X, NumericVector g, NumericVector y, NumericVector xx, NumericVector E, NumericVector L, int N, double Ve);
RcppExport SEXP NAM_SAMP2(SEXP XSEXP, SEXP gSEXP, SEXP ySEXP, SEXP xxSEXP, SEXP ESEXP, SEXP LSEXP, SEXP NSEXP, SEXP VeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g(gSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type Ve(VeSEXP);
    SAMP2(X, g, y, xx, E, L, N, Ve);
    return R_NilValue;
END_RCPP
}
// timesMatrix
NumericMatrix timesMatrix(NumericMatrix ma1, NumericVector h, NumericMatrix ma2, int rows, int cols);
RcppExport SEXP NAM_timesMatrix(SEXP ma1SEXP, SEXP hSEXP, SEXP ma2SEXP, SEXP rowsSEXP, SEXP colsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ma1(ma1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ma2(ma2SEXP);
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    rcpp_result_gen = Rcpp::wrap(timesMatrix(ma1, h, ma2, rows, cols));
    return rcpp_result_gen;
END_RCPP
}
// timesVec
NumericMatrix timesVec(NumericVector aa, NumericVector h, NumericMatrix bb, int q);
RcppExport SEXP NAM_timesVec(SEXP aaSEXP, SEXP hSEXP, SEXP bbSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type aa(aaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(timesVec(aa, h, bb, q));
    return rcpp_result_gen;
END_RCPP
}
