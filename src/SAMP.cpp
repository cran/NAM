#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void SAMP(NumericMatrix C, NumericVector g, NumericVector r, int N, double Ve)
{
  RNGScope scope;
  for(int i=0; i<N; i++){        
    g[i]=R::rnorm(((r[i]-sum(C(i,_)*g)+C(i,i)*g(i))/C(i,i)),sqrt(Ve/C(i,i)));
  } 
}
