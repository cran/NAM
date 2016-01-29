#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector funX(NumericVector col, int finsiz) 
{
  int scol = col.size();
  int i;
  int ct = 0;
  NumericVector vec2(finsiz);
  
  for(i=0; i < scol; i++)
  {
    if(col[i] == 0 || col[i] == 2)
    {
      vec2[ct] = 2;
    }
    else if(col[i] == 1)
    {
      vec2[ct] = 1;
      ct++;
      vec2[ct] = 1;
    }
    ct++;    
  }
  return (vec2);
}
