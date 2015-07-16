#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
void SAMP2(NumericMatrix X, NumericVector g, NumericVector y, NumericVector xx, NumericVector E, NumericVector L, int N, double Ve)
{
  RNGScope scope; 
  double G;
  for(int i=0; i<N; i++){
    G = g[i];
    g[i] = R::rnorm(
    (sum(X(_,i)*E) + xx(i)*G)/(xx(i)+L(i)),
    sqrt(Ve/(xx(i)+L(i))));
    E = E - X(_,i) * (g[i]-G);
  }
}
