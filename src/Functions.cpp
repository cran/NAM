#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP emBA(NumericVector y, NumericMatrix gen, double df = 4, double R2 = 0.5){
  int it = 200;
  int p = gen.ncol();
  int n = gen.nrow();
  double ve = 1;
  double va = 1;
  double vy = var(y);
  NumericVector xx(p);
  NumericVector vx(p);
  for(int k=0; k<p; k++){
    xx[k] = sum(gen(_,k)*gen(_,k));
    vx[k] = var(gen(_,k));
  }
  NumericVector b(p);
  NumericVector vb = b+va;
  NumericVector Lmb = ve/vb;
  double MSx = sum(vx);
  double Sb = R2*(df+2)*vy/MSx;
  double Se = (1-R2)*(df+2)*vy;
  double shape_hp  = 1.1;
  double rate_hp = (shape_hp-1)/Sb;
  double mu = mean(y);
  NumericVector e = y-mu;
  double b0,eM,h2,Sb_j;
  for(int i=0; i<it; i++){
    Sb_j = (p*df/2+shape_hp)/(sum(1/vb)/2+rate_hp);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b[j] = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]);
      e = e-gen(_,j)*(b[j]-b0);
      vb[j] = (Sb_j+b[j]*b[j]+va)/(df+2)+ve/(xx[j]+Lmb[j]);
    }
    va = sum(b*b)/p;
    ve = (sum(e*e)+Se)/(n+df);
    Lmb = ve/vb;
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
  }
  double vg = sum(vb);
  h2 = vg/(vg+ve);
  NumericVector fit(n);
  for(int k=0; k<n; k++){
    fit[k] = sum(gen(k,_)*b)+mu;
  }
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("Va") = vg,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emBB(NumericVector y, NumericMatrix gen, double df = 4, double R2 = 0.5, double Pi = 0.7){
  R2 = R2/5;
  int it = 200;
  df = df*10;
  int p = gen.ncol();
  int n = gen.nrow();
  double va = 1;
  double ve = 1;
  NumericVector d(p);
  NumericVector b(p);
  NumericVector vb = b+va;
  NumericVector Lmb = ve/vb;
  double vy = var(y);
  if(Pi>0.5){
    Pi = 1-Pi;
  } 
  NumericVector xx(p);
  NumericVector vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(gen(_,i)*gen(_,i));
    vx[i] = var(gen(_,i));
  }
  double MSx = sum(vx)*Pi;
  double Sb = R2*(df+2)*vy/MSx;
  double Se = (1-R2)*(df+2)*vy;
  double shape_hp  = 1.1;
  double rate_hp = (shape_hp-1)/Sb;
  double mu = mean(y);
  NumericVector e = y-mu;
  NumericVector e1(n);
  NumericVector e2(n);
  double b0,b1,LR,eM,h2,Sb_j,C;
  double Pi0 = (1-Pi)/Pi;
  double MD = Pi;
  for(int i=0; i<it; i++){
    C = -0.5/ve;
    Sb_j = (p*df/2+shape_hp)/(sum(1/vb)/2+rate_hp);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]);
      e1 = e-gen(_,j)*(b1-b0);
      e2 = e-gen(_,j)*(0-b0);
      LR = Pi0*exp(C*(sum(e2*e2)-sum(e1*e1)));
      d[j] = (1/(1+LR));
      b[j] = b1*d[j]/MD;
      vb[j] = (Sb_j+b[j]*b[j]+va)/(df+2);
      e = e - gen(_,j)*(b1-b0);
    }
    MD = max(d);
    va = sum(b*b)/p;
    ve = (sum(e*e)+Se)/(n+df);
    Lmb = ve/vb;
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
  }
  double vg = sum(vb)/Pi;
  h2 = vg/(vg+ve);
  NumericVector fit(n);
  for(int k=0; k<n; k++){
    fit[k] = sum(gen(k,_)*b)+mu;
  }
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("Va") = vg,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emBC(NumericVector y, NumericMatrix gen, double df = 4, double R2 = 0.5, double Pi = 0.7){
  R2 = R2/5;
  int it = 200;
  int p = gen.ncol();
  int n = gen.nrow();
  double Lmb = 1;
  double vb = 1;
  double va = 1;
  double ve = 1;
  double vy = var(y);
  if(Pi>0.5){
    Pi = 1-Pi;
  } 
  NumericVector xx(p);
  NumericVector vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(gen(_,i)*gen(_,i));
    vx[i] = var(gen(_,i));
  }
  double MSx = sum(vx)*Pi;
  double Sb = R2*(df+2)*vy/MSx;
  double Se = (1-R2)*(df+2)*vy;
  double mu = mean(y);
  NumericVector b(p);
  NumericVector d(p);
  NumericVector e = y-mu;
  NumericVector e1(n);
  NumericVector e2(n);
  double b0,b1,LR,eM,h2,C;
  double Pi0 = (1-Pi)/Pi;
  double MD = Pi;
  for(int i=0; i<it; i++){
    C = -0.5/ve;
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);
      e1 = e-gen(_,j)*(b1-b0);
      e2 = e-gen(_,j)*(0-b0);
      LR = Pi0*exp(C*(sum(e2*e2)-sum(e1*e1)));
      d[j] = (1/(1+LR));
      b[j] = b1*d[j];
      e = e - gen(_,j)*(b1-b0);
    }
    MD = max(d);
    b = b/MD;
    vb = (sum(b*b)+Sb)/(p+df);
    ve = (sum(e*e)+Se)/(n+df);
    Lmb = ve/vb;
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
  }
  va = mean(xx)*(sum(b*b)+Sb)/(p+df)/Pi;
  h2 = va/(va+ve);
  NumericVector fit(n);
  for(int k=0; k<n; k++){
    fit[k] = sum(gen(k,_)*b)+mu;
  }
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("Va") = va,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emRR(NumericVector y, NumericMatrix gen, double df = 4, double R2 = 0.5){
  int it = 200;
  df = df*10000;
  int p = gen.ncol();
  int n = gen.nrow();
  double Lmb = 1;
  double vb = 1;
  double va = 1;
  double ve = 1;
  double vy = var(y);
  NumericVector xx(p);
  NumericVector vx(p);
  for(int k=0; k<p; k++){
    xx[k] = sum(gen(_,k)*gen(_,k));
    vx[k] = var(gen(_,k));
  }
  double MSx = sum(vx);
  double Sb = R2*(df+2)*vy/MSx;
  double Se = (1-R2)*(df+2)*vy;
  double mu = mean(y);
  NumericVector b(p);
  NumericVector e = y-mu;
  double b0,eM,h2;
  for(int i=0; i<it; i++){
    for(int j=0; j<p; j++){
      b0 = b[j];
      b[j] = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);
      e = e-gen(_,j)*(b[j]-b0);
    }
    vb = (sum(b*b)+Sb)/(p+df+2);
    ve = (sum(e*e)+Se)/(n+df+2);
    Lmb = ve/vb;
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
  }
  va = sqrt(vb)*mean(xx);
  ve = sqrt(ve);
  h2 = (va/(va+ve))*(va/(va+ve));
  NumericVector fit(n);
  for(int k=0; k<n; k++){
    fit[k] = sum(gen(k,_)*b)+mu;
  }
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("Va") = va,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emBL(NumericVector y, NumericMatrix gen, double R2 = 0.5, double alpha = 0.02){
  int it = 200;
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
  h2 = 1 - var(e)/var(y);
  // Output
  return List::create(Named("mu") = mu, Named("b") = b, Named("hat") = fit, Named("h2") = h2); }

// [[Rcpp::export]]
SEXP emDE(NumericVector y, NumericMatrix gen, double R2 = 0.5){
  // Convergence criteria
  int maxit = 300;
  double tol = 10e-6;
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  double b0,eM;
  double mu = mean(y);
  NumericVector e = y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
  NumericVector vx(p);
  for(int k=0; k<p; k++){vx[k] = var(gen(_,k));}
  double cxx = sum(vx)*(1-R2)/R2;
  // Regulation coefficients
  double Ve;
  NumericVector Vb(p);
  NumericVector b(p);
  NumericVector Lmb = p+cxx;
  double b1;
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      // Ordinary Least Square
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(Lmb(j)+xx(j));
      b[j] = b1;
      // Residuals update
      e = e-gen(_,j)*(b1-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Variance components
    Ve = sum(e*y)/(n-1);
    Vb = b*b+(Ve/(xx+Lmb));
    Lmb = sqrt(cxx*Ve/Vb);
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}
  }
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  // Output
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("Vb")=Vb,
                      Named("Ve")=Ve,
                      Named("h2")=sum(Vb)/(sum(Vb)+Ve));
}

// [[Rcpp::export]]
SEXP emEN(NumericVector y, NumericMatrix gen, double R2 = 0.5, double alpha = 0.02){
  // Convergence criteria
  int maxit = 300;
  double tol = 10e-11;
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  NumericVector b(p);
  double b0,eM;
  double mu = mean(y);
  NumericVector e = y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
  NumericVector vx(p);
  for(int k=0; k<p; k++){vx[k] = var(gen(_,k));}
  double cxx = sum(vx)*(1-R2)/R2;
  // Regulation coefficients
  double Ve, Va;
  double Sy = sd(y);
  double Lmb = cxx;
  double Lmb1 = 0.5*Lmb*alpha*Sy;
  double Lmb2 = Lmb*(1-alpha);
  double trAC22 = sum(1/(xx+Lmb));
  double OLS, b1;
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      // Ordinary Least Square
      b0 = b[j];
      OLS = (sum(gen(_,j)*e)+xx[j]*b0);
      // Regularization of positive OLS
      if(OLS>0){
        // Regularization of positive OLS
        b1 = (OLS-Lmb1)/(Lmb2+xx(j));if(b1<0){b1=0;};
      }else{
        // Regularization of negative OLS
        b1 = (OLS+Lmb1)/(Lmb2+xx(j));if(b1>0){b1=0;};
      }
      b[j] = b1;
      // Residuals update
      e = e-gen(_,j)*(b1-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Variance components
    Ve = sum(e*y)/(n-1);
    Va = (sum(b*b)+trAC22*Ve)/p;
    Lmb = Ve/Va;
    Lmb1 = 0.5*Lmb*alpha*Sy;
    Lmb2 = Lmb*(1-alpha);
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}
  }
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  // Output
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("Va")=Va*cxx,
                      Named("Ve")=Ve,
                      Named("h2")=Va*cxx/(Va*cxx+Ve));
}

// [[Rcpp::export]]
int calcSize(NumericVector col, NumericVector fam){
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

// [[Rcpp::export]]
NumericVector funI(NumericVector col, int fam, int finsiz, int f){
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

// [[Rcpp::export]]
NumericVector funX(NumericVector col, int finsiz){
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

// [[Rcpp::export]]
void gs(NumericMatrix C, NumericVector g, NumericVector r, int N){
  for(int i=0; i<N; i++){        
    g[i]=((r[i]-sum(C(i,_)*g)+C(i,i)*g(i))/C(i,i));
  } 
}

// [[Rcpp::export]]
NumericVector inputRow(NumericVector x, int exp, int n){  
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

// [[Rcpp::export]]
SEXP KMUP(NumericMatrix X, NumericVector b,
          NumericVector d, NumericVector xx,
          NumericVector E, NumericVector L,
          double Ve, double pi){
  
  RNGScope scope;
  int p = X.ncol();
  NumericVector g = b;
  NumericVector e1 = E;
  NumericVector e2 = E;
  double G,G0,Gp,D,pj,LR;
  double Cons = -0.5/Ve;
  double Pi0 = pi/(1-pi);
  
  for(int i=0; i<p; i++){
    G0 = g[i];
    G = R::rnorm((sum(X(_,i)*E)+G0)/(xx(i)+L(i)),
                 sqrt(Ve/(xx(i)+L(i))));
    Gp = G*pi;
    e1 = E - X(_,i)*(G-G0);
    if(pi>0){
      e2 = E - X(_,i)*(Gp-G0);
      LR = Pi0*exp(Cons*(sum(e1*e1)-sum(e2*e2)));
      pj = 1-1/(1+LR);
      D = R::rbinom(1,pj);
      if(D==0){
        d[i] = 0;
        g[i] = Gp;
        E = e2;
      }else{
        d[i] = 1;
        g[i] = G;
        E = e1;
      }
    }else{
      d[i] = 1;
      g[i] = G;
      E = e1;
    }
  }
  return List::create(Named("b") = g,
                      Named("d") = d,
                      Named("e") = E);
}

// [[Rcpp::export]]
SEXP KMUP2(NumericMatrix X, NumericVector Use, NumericVector b,
           NumericVector d, NumericVector xx, NumericVector E,
           NumericVector L, double Ve, double pi){
  
  RNGScope scope;
  int p = X.ncol();
  int n0 = X.nrow();
  int n = Use.size();
  double bg = n0/n;
  NumericVector E0(n);
  NumericVector H(n);
  for(int k=0; k<n; k++){
    E0[k] = E[Use[k]];
  }
  NumericVector g = b;
  NumericVector e1 = E0;
  NumericVector e2 = E0;
  double G,G0,Gp,D,pj,LR;
  double Cons = -0.5/Ve;
  double Pi0 = pi/(1-pi);
  
  for(int i=0; i<p; i++){
    G0 = g[i];
    for(int x=0; x<n; x++){
      H[x] = X(Use[x],i);
    }
    G = R::rnorm((sum(H*E0)+G0)/(xx(i)*bg+L(i)),
                 sqrt(Ve/(xx(i)*bg+L(i))));
    Gp = G*pi;
    e1 = E0 - H*(G-G0);
    if(pi>0){
      e2 = E0 - H*(Gp-G0);
      LR = Pi0*exp(Cons*(sum(e1*e1)-sum(e2*e2)));
      pj = 1-1/(1+LR);
      D = R::rbinom(1,pj);
      if(D==0){
        d[i] = 0;
        g[i] = Gp;
        E0 = e2;
      }else{
        d[i] = 1;
        g[i] = G;
        E0 = e1;
      }
    }else{
      d[i] = 1;
      g[i] = G;
      E0 = e1;
    }
  }
  
  
  return List::create(Named("b") = g,
                      Named("d") = d,
                      Named("e") = E0);
}

// [[Rcpp::export]]
void SAMP(NumericMatrix C, NumericVector g, NumericVector r, int N, double Ve){
  RNGScope scope;
  for(int i=0; i<N; i++){        
    g[i]=R::rnorm(((r[i]-sum(C(i,_)*g)+C(i,i)*g(i))/C(i,i)),sqrt(Ve/C(i,i)));
  } 
}

// [[Rcpp::export]]
void SAMP2(NumericMatrix X, NumericVector g, NumericVector y, NumericVector xx, NumericVector E, NumericVector L, int N, double Ve){
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

// [[Rcpp::export]]
NumericMatrix timesMatrix(NumericMatrix ma1, NumericVector h, NumericMatrix ma2, int rows, int cols){
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

// [[Rcpp::export]]
NumericMatrix timesVec(NumericVector aa, NumericVector h, NumericMatrix bb, int q){
  int i;
  NumericMatrix resul(q,1);
  
  for(i = 0; i < q; i++)
  {
    resul[i] = sum(aa*h*bb(_,i));
  }
  
  return(resul);
}

// [[Rcpp::export]]
void CNT(NumericMatrix X){for(int j=0;j<X.ncol();j++){X(_,j)=X(_,j)-mean(X(_,j));}}

// [[Rcpp::export]]
SEXP MSX(NumericMatrix X){
  int p = X.ncol(); int n = X.nrow(); double m; NumericVector xx(p); NumericVector sx(p);
  for(int k=0; k<p; k++){ xx[k] = sum(X(_,k)*X(_,k)); m = sum(X(_,k)); sx[k] = m*m/n; }
  double cxx = sum(xx-sx)/(n-1); return List::create(Named("cxx")=cxx,Named("xx")=xx);}

// [[Rcpp::export]]
void IMP(NumericMatrix X){;int p = X.ncol(); int n = X.nrow();
LogicalVector MIS(n); NumericVector x(n); NumericVector z; double EXP;
for(int j=0; j<p; j++){;if(is_true(any(is_na(X(_,j))))){
  x = X(_,j); MIS = is_na(x);z = x[!MIS]; EXP = mean(z);
  X(_,j) = ifelse(MIS,EXP,x);};};};

// [[Rcpp::export]]
SEXP NOR(NumericVector y, NumericMatrix X, double cxx, NumericVector xx, int maxit = 50, double tol = 10e-6){
  int p = X.ncol(); int n = X.nrow(); double Ve,b0,b1,eM; double mu = mean(y); double Va=1;
  NumericVector b(p), bc(p), fit(n), E(n); NumericVector e = y-mu; double Lmb = cxx;
  int numit = 0; double cnv = 1; while(numit<maxit){bc = b+0; for(int j=0; j<p; j++){
    b0 = b[j]; b1 = (sum(X(_,j)*e)+xx[j]*b0)/(Lmb+xx(j)); b[j]=b1; e = e-X(_,j)*(b1-b0);}
  eM = mean(e); mu = mu+eM; e = e-eM; Ve = sum(e*y)/(n-1); Va = var(b)+mean(Ve/(xx+Lmb));
  Lmb = sqrt(cxx*Ve/Va); ++numit; cnv = sum(abs(bc-b)); if( cnv<tol ){break;};}
  for(int k=0; k<n; k++){ fit[k] = sum(X(k,_)*b); E[k] = y[k]-fit[k];}; Va=Va*cxx;
  return List::create(Named("b")=b,Named("v")=Va,Named("h")=fit,Named("e")=E);}

// [[Rcpp::export]]
NumericMatrix GAU(NumericMatrix X){
  int n = X.nrow(); NumericVector D; NumericMatrix K(n,n); double d2, md;
  for(int i=0; i<n; i++){; for(int j=0; j<n; j++){
    if(i==j){ K(i,j)=0; }else if(j>i){; D = X(i,_)-X(j,_);
    d2 = sum(D*D); d2 = d2*d2; K(i,j)=d2; K(j,i)=d2; }}}; md = mean(K);
    for(int i=0; i<n; i++){K(i,_) = exp(-K(i,_)/md);} return K;}

// [[Rcpp::export]]
NumericVector SPC(NumericVector y, NumericVector blk, NumericVector row, NumericVector col, int rN=3, int cN=1){
  int n = y.size(); NumericVector Cov(n), Phe(n), Obs(n);
  for(int i=0; i<n; i++){; for(int j=0; j<n; j++){
    if( (i!=j) & (blk[i]==blk[j]) & (abs(row[i]-row[j])<=rN) & (abs(col[i]-col[j])<=cN) ){
      Phe[j] = y[j]; Obs[j] = 1; }else{ Phe[j] = 0; Obs[j] = 0; }}
  Cov[i] = sum(Phe)/sum(Obs);} return Cov;}