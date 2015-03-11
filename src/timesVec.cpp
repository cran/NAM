#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix timesVec(NumericVector aa, NumericVector h, NumericMatrix bb, int q)
{
   int i;
   NumericMatrix resul(q,1);
   
   for(i = 0; i < q; i++)
   {
     resul[i] = sum(aa*h*bb(_,i));
   }
   
   return(resul);
}
