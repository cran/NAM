#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void gs(NumericMatrix C, NumericVector g, NumericVector r, int N)
{
  for(int i=0; i<N; i++){        
    g[i]=((r[i]-sum(C(i,_)*g)+C(i,i)*g(i))/C(i,i));
  } 
}
