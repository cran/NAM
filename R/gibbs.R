
gibbs = function(y,X=NULL,Z=NULL,iK=NULL,Iter=1500,Burn=500,Thin=4,DF=5){
    
  # Default for X; changing X to matrix if it is a formulas
  if(is.null(X)) X=matrix(1,length(y),1)
  if(class(X)=="formula"){
    X=model.frame(X)
    Fixes=ncol(X)
    XX=matrix(1,length(y),1)
    for(var in 1:Fixes) XX=cbind(XX,model.matrix(~X[,var]-1))
    X=XX
    rm(XX,Fixes)
  }
  
  # Defaults for Z, and making "Z formula" as "Z list"
  if(class(Z)=="formula") {
    Z=model.frame(Z)
    Randoms=ncol(Z)
    ZZ=list()
    for(var in 1:Randoms) ZZ[[var]]=model.matrix(~Z[,var]-1)
    Z=ZZ
    rm(ZZ,Randoms)
  }
  if(is.null(Z)) Z=list(diag(length(y)))
  
  # Defaults for null and incomplete iK
  if(is.null(iK)){
    iK=list()
    Randoms=length(Z)
    for(var in 1:Randoms) iK[[var]]=diag(ncol(Z[[var]]))
  }
  if(class(iK)=="matrix") iK=list(iK)
  if(length(Z)!=length(iK)){
    a=length(Z)
    b=length(iK)
    if(a>b) for(K in 1:(a-b)) iK[[(K+b)]]=diag(ncol(Z[[K]]))
    if(b>a) for(K in 1:(b-a)) Z[[(K+a)]]=diag(ncol(iK[[K]]))
    rm(a,b)
  }
  
  # Predictiors should not have missing values
  if(any(is.na(X))|any(is.na(unlist(Z)))) stop("Predictors with missing values not allowed")
  
  # Add logit link function
  # Add residual variation
  
  # Gibb Sampling 
  # Garcia-Cortes, L. A. & Sorensen, D. (1996). GSE, 28(1), 121-126.
  # Sorensen, D., & Gianola, D. (2002). Springer.
  # Prior solution: de los Campos et al (2013). Genetics, 193(2), 327-345.
  
  # Thinning - which Markov Chains are going to be stored
  THIN = seq(Burn,Iter,Thin)
  
  # Some parameters with similar notation from the book
  nx = ncol(X)
  Randoms = length(Z) # number of random variables
  q = rep(0,Randoms); for(i in 1:Randoms) q[i]=ncol(Z[[i]])  
  N = nx+sum(q)
  
  # Qs1 and Qs2 regard where each random variable starts and ends, respectively
  Qs0 = c(nx,q)
  Variables = length(Qs0)
  Qs1 = Qs2 = rep(0,Variables)
  for(i in 1:Variables){
    Qs1[i]=max(Qs2)+1
    Qs2[i]=Qs1[i]+Qs0[i]-1
  }
  # Starting values for the variance components
  Ve = 1
  Va = lambda = rep(1,Randoms)
  
  # Linear system described as: WW+Sigma = Cg = r
  W = X
  for(i in 1:Randoms) W=cbind(W,Z[[i]])
  # MISSING
  W1=W
  if(any(is.na(y))){MIS = which(is.na(y));W=W[-MIS,];y=y[-MIS]} 
  n = length(y)
  # Keeping on
  r = crossprod(W,y) # Need adjust to accept W.iR.y
  WW = (crossprod(W)) # Need adjust to accept W.iR.W
  # Covariance Matrix
  Sigma = matrix(0,N,N)
  for(i in 1:Randoms) Sigma[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]] = iK[[i]]*lambda[i]
  # Matching WW and Sigma
  C = WW+Sigma
  g = rep(0,N)
  
  # Saving space for the posterior
  include = 0
  POSTg = matrix(0,ncol = length(THIN), N)
  POSTv = matrix(0,ncol = length(THIN), nrow = (Randoms+1))
  
  # Prior degrees of freedom Prior solution for shape/variances
  df0 = DF
  S0=rep(0,Randoms) # S0 has analystical solution (De los Campos et al 2013)
  for(i in 1:Randoms) 
    S0[i]=(var(y,na.rm=T)*0.5)/mean(
      (t((t(C[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]])-
            colMeans(C[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]]))))^2 )
  
  # Saving memory for some vectors
  e = rep(0,N)
  
  # Progression Bar
  pb=txtProgressBar(style=3)
  
  # LOOP
  for(iteration in 1:Iter){
    
    # Sampling hyper-priors
    S0a = runif(Randoms,S0*0.66,S0*1.33)
    df0a = runif(Randoms,df0*0.66,df0*1.33)
    dfu = q + df0a
    S0b = runif(1,0.5,1.2)
    df0b = runif(1,df0*0.66,df0*1.33)
    dfe = n + df0b
    
    # Random variance
    for(i in 1:Randoms)
      # (ZiAZ+S0v0)/x2(v)
      Va[i] = (sum(crossprod(g[Qs1[i+1]:Qs2[i+1]],iK[[i]])*(g[Qs1[i+1]:Qs2[i+1]]))+  
                 S0a[i]*df0a[i])/rchisq(1,df=dfu[i])
    
    # Residual variance
    e = y - W%*%g
    Ve = (crossprod(e)+S0b*df0b) / rchisq(1,df=dfe)
    
    # Ve/Va
    lambda = Ve/Va
    
    # Updating C
    for(i in 1:Randoms)
      Sigma[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]] = iK[[i]]*lambda[i]
    C = WW+Sigma
    
    SAMP(C,g,r,N,Ve) # the C++ SAMP updates "g" and doesn't return anything
    
    if(is.element(iteration,THIN)){
      include = include + 1;
      POSTg[,include] = g;
      POSTv[1:Randoms,include] = Va;
      POSTv[(Randoms+1),include] = Ve;
    }  
    # Advance progression bar
    setTxtProgressBar(pb,iteration/Iter)
  }
  # End progression bar
  close(pb)
  
  # Mode Function - Venter J.H. (1967). Ann. Math. Statist., 38(5):1446-1455.
  moda=function (x){
    it=5;ny=length(x);k=ceiling(ny/2)-1; while(it>1){
      y=sort(x); inf=y[1:(ny-k)]; sup=y[(k+1):ny]
      diffs=sup-inf; i=min(which(diffs==min(diffs)))
      M=median(y[i:(i+k)]); it=it-1}; return(M)}
  
  rownames(POSTg)=paste("b",0:(N-1),sep="")
  for(i in 1:Randoms) rownames(POSTg)[Qs1[i+1]:Qs2[i+1]] = paste("u",i,".",1:Qs0[i+1],sep="")
  
  # Mean and Mode Posterior
  Post.B = apply(POSTg,1,moda)
  Mean.B = apply(POSTg,1,mean)
  Post.VC = c(apply(POSTv,1,moda))
  names(Post.VC)=c(paste("Va",1:Randoms,sep=""),"Ve")
  
  # List of Coefficients
  Coefficients = list()
  Coefficients[[1]] = POSTg[1:nx,]
  for(i in 1:Randoms) Coefficients[[i+1]]=POSTg[(Qs1[i+1]:Qs2[i+1]),]
  names(Coefficients)[1]="Fixed"
  for(i in 1:Randoms) names(Coefficients)[i+1]=paste("Random",i,sep="")
  
  # Return
  return( list("Coef.mode" = Post.B,
               "Coef.mean" = Mean.B,
               "VC.mode" = Post.VC,
               "Posterior.Coef" = Coefficients,
               "Posterior.VC" = POSTv,
               "Fit.mean" = W1%*%Post.B,
               "Fit.mode" = W1%*%Mean.B
  ) )
  
}
