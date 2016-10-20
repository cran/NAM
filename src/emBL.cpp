#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP emBL(NumericVector y, NumericMatrix gen, int it = 200, double R2 = 0.5, double alpha = 0.025){
  double h2 = R2;  
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  NumericVector b(p);
  double b0,eM,Half_L2;
  double mu = mean(y);
  NumericVector e = y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
  // Regulation coefficients
  double cxx = mean(xx);
  double Lmb1 = cxx*((1-h2)/h2)*alpha*0.5;
  double Lmb2 = cxx*((1-h2)/h2)*(1-alpha);
  double OLS, G;
  // Loop
  for(int i=0; i<it; i++){
    // Regression coefficients loop
    for(int j=0; j<p; j++){
    // Ordinary Least Square
    b0 = b[j];
    OLS = (sum(gen(_,j)*e)+xx[j]*b0);
    Half_L2 = 0.5*OLS/(xx[j]+cxx);
    // Regularization of positive OLS
    if(OLS>0){
      G = 0.5*(OLS-Lmb1)/(Lmb2+xx(j));
      if(G>0){b[j] = G+Half_L2;}else{b[j] = Half_L2;}
    }else{
    // Regularization of negative OLS
      G = 0.5*(OLS+Lmb1)/(Lmb2+xx(j));
      if(G<0){b[j] = G+Half_L2;}else{b[j] = Half_L2;}
    }
    // Residuals update
    e = e-gen(_,j)*(b[j]-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;}
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  // Output
  return List::create(Named("mu") = mu, Named("b") = b, Named("hat") = fit); }
