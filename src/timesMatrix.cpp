#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix timesMatrix(NumericMatrix ma1, NumericVector h, NumericMatrix ma2, int rows, int cols)
{
  int i;
  int j;
  NumericMatrix resul(rows,cols);
  
  for(i = 0; i < rows; i++)
  {
    for(j = 0; j < cols; j++)
    {
      resul(i,j) = sum(ma1(_,i)*h*ma2(_,j)); 
    }
  }
  
   return (resul);
}
