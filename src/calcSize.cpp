#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int calcSize(NumericVector col, NumericVector fam) 
{
  int scol = col.size();
  int i;
  int finsiz = 0;
  
  for(i=0; i < scol; i++)
  {
    if(col[i] == 1)
    {
      finsiz = finsiz + 2;
    }
    else
    {
      finsiz = finsiz + 1;
    }
  }
  return (finsiz);
}
