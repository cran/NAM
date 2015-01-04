#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector funI(NumericVector col, int fam, int finsiz, int f) 
{
  int scol = col.size();
  NumericVector vec(finsiz);
  int ct = 0;
  //NumericVector vec(finsiz);
  
  ct = 0;
  
  for(int i=0; i < scol; i++)
  {
    if(col[i] == 2)
    {
      vec[ct] = (i) * f + 1;
    }
    else if(col[i] == 1)
    {
      vec[ct] = (i) * f + 1;
    //  printf("Value I: %d\n", vec[ct]);
      ct++;
      vec[ct] = (i) * f + 1 + fam;
    }
    else if(col[i] == 0)
    {
      vec[ct] = (i) * f + 1 + fam;
    }
   // printf("Value I: %d\n", vec[ct]);
    ct++;
  }
  
   return vec;
}
