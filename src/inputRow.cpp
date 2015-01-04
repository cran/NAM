#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector inputRow(NumericVector x, int exp, int n) 
{  
  int i;
  
   // IF THE FIRST IS MISSING
   if(x[0] == 5 && x[1] == 5)
   {
     x[0] = exp;
   }
   else if(x[0] == 5)
   {
     x[0] = x[1];
   }
   // IF ANY OTHER IS MISSING
   for(i = 1; i < n; i++)
   {
     if(x[i] == 5 && x[i+1] == 5)
     {
       x[i] = x[i-1];
     }
     else if(x[i] == 5 && (x[i-1] == x[i+1]))
     {
       x[i] = x[i-1];
     }
     else if(x[i] == 5 && x[i-1] != x[i+1])
     {
       x[i] = 1;
     }
   
   }
       
   return x;
}
