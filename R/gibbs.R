
gibbs = function(y,Z=NULL,X=NULL,iK=NULL,iR=NULL,Iter=1500,Burn=500,Thin=4,DF=5,S=1,GSRU=FALSE){
  anyNA = function(x) any(is.na(x))  
  # Default for X; changing X to matrix if it is a formulas
  VY = var(y,na.rm=T)
  if(is.null(X)) X=matrix(1,length(y),1)
  if(class(X)=="formula"){
    X=model.frame(X)
    Fixes=ncol(X)
    XX=matrix(1,length(y),1)
    for(var in 1:Fixes) XX=cbind(XX,model.matrix(~X[,var]-1))
    X=XX
    rm(XX,Fixes)
  }
  
  # Defaults of Z: making "NULL","formula" and "matrix" as "list"
  if(is.null(Z)&is.null(iK)) stop("Either Z or iK must be specified")
  if(is.null(Z)) Z=list(diag(length(y)))
  if(class(Z)=="matrix") Z = list(Z)
  if(class(Z)=="formula") {
    Z=model.frame(Z)
    Randoms=ncol(Z)
    ZZ=list()
    for(var in 1:Randoms) ZZ[[var]]=model.matrix(~Z[,var]-1)
    Z=ZZ
    rm(ZZ,Randoms)
  }

  if(!GSRU){
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
  
  # Some parameters with notation from the book
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
  if(any(is.na(y))){
    MIS = which(is.na(y))
    W=W[-MIS,]
    y=y[-MIS]
    if(!is.null(iR)) iR=iR[-MIS,-MIS]
  }
  n = length(y)
  # Keeping on
  # GSRU does not deal with WW or iK
  if(!GSRU){
  if(is.null(iR)){
    r = crossprod(W,y)
    WW = (crossprod(W))
  }else{
    r = crossprod(W,iR)%*%y
    WW = (crossprod(W,iR))%*%W
  }
  # Covariance Matrix
  Sigma = matrix(0,N,N)
  for(i in 1:Randoms) Sigma[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]] = iK[[i]]*lambda[i]
  # Matching WW and Sigma
  C = WW+Sigma
  }else{
    L = rep(0,ncol(W))
    for(i in 1:Randoms) L[Qs1[i+1]:Qs2[i+1]] = lambda[i]
    xx = colSums(W^2)
  }
  g = rep(0,N)
  # Saving space for the posterior
  include = 0
  POSTg = matrix(0,ncol = length(THIN), N)
  POSTv = matrix(0,ncol = length(THIN), nrow = (Randoms+1))
  
  # Hperpriors: Degrees of freedom (DF) and Shape (S)
  df0 = DF
  
  if(is.null(S)){
    S0=rep(0,Randoms) # S0 has analystical solution (De los Campos et al 2013)
    if(!GSRU){
    for(i in 1:Randoms) 
      S0[i]=(var(y,na.rm=T)*0.5)/mean(
        (t((t(C[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]])-
              colMeans(C[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]]))))^2 )
    }else{
      for(i in 1:Randoms) 
        S0[i]=(var(y,na.rm=T)*0.5)/
        mean( (t((t(W[,Qs1[i+1]:Qs2[i+1]])-colMeans(W[,Qs1[i+1]:Qs2[i+1]]))))^2 )
    }
    
  }else{S0=rep(S,Randoms)}
    
  # Saving memory for some vectors
  e = rep(0,N)
  
  # Progression Bar
  pb=txtProgressBar(style=3)
  
  # LOOP
  for(iteration in 1:Iter){
    
    # Sampling hyper-priors
    S0a = runif(Randoms,S0*0.5,S0*1.5)
    df0a = runif(Randoms,df0*0.5,df0*1.5)
    dfu = q + df0a
    if(is.null(S)){ S0b = runif(1,0.0001,5) }else{ S0b = runif(1,S*0.5,S*1.5) }
    df0b = runif(1,min(2,df0)*0.5,min(2,df0)*1.5)
    dfe = n + df0b
    
    # Random variance
    for(i in 1:Randoms){
      # (ZiAZ+S0v0)/x2(v)
      if(!GSRU){
      Va[i] = (sum(crossprod(g[Qs1[i+1]:Qs2[i+1]],iK[[i]])*(g[Qs1[i+1]:Qs2[i+1]]))+  
                 S0a[i]*df0a[i])/rchisq(1,df=dfu[i])
      }else{
        Va[i] = (sum(crossprod(g[Qs1[i+1]:Qs2[i+1]]))+S0a[i]*df0a[i])/rchisq(1,df=dfu[i])
      }
    }
      
    # Residual variance
    e = y - W%*%g
    Ve = (crossprod(e)+S0b*df0b) / rchisq(1,df=dfe)
    
    # Ve/Va
    lambda = Ve/Va
    
    # Updating C
    if(!GSRU){
    for(i in 1:Randoms) Sigma[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]] = iK[[i]]*lambda[i]
    C = WW+Sigma
    }
    
    # the C++ SAMP updates "g" and doesn't return anything
    if(!GSRU){
      SAMP(C,g,r,N,Ve)
    } else {
      SAMP2(W,g,y,xx,e,L,N,Ve) 
    }
        
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
  Mean.B = apply(POSTg,1,mean)
  Post.VC = c(apply(POSTv,1,moda))
  names(Post.VC) = c(paste("Va",1:Randoms,sep=""),"Ve")
  rownames(POSTv) = c(paste("Va",1:Randoms,sep=""),"Ve")
  
  # List of Coefficients
  Coefficients = list()
  Coefficients[[1]] = POSTg[1:nx,]
  for(i in 1:Randoms) Coefficients[[i+1]]=POSTg[(Qs1[i+1]:Qs2[i+1]),]
  names(Coefficients)[1]="Fixed"
  for(i in 1:Randoms) names(Coefficients)[i+1]=paste("Random",i,sep="")
  
  
  RESULTS = list(
                 "Coef.estimate" = Mean.B,
                 "VC.estimate" = Post.VC,
                 "Posterior.Coef" = Coefficients,
                 "Posterior.VC" = POSTv,
                 "Fit.mean" = W1%*%Mean.B
                 )
  
  class(RESULTS) = "gibbs"
  
  # Return
  return( RESULTS )
  
}

ml = function(y,Z=NULL,X=NULL,iK=NULL,iR=NULL,DF=-2,S=0){
  
  anyNA = function(x) any(is.na(x))
  
  # Default for X; changing X to matrix if it is a formulas
  VY = var(y,na.rm=T)
  if(is.null(X)) X=matrix(1,length(y),1)
  if(class(X)=="formula"){
    X=model.frame(X)
    Fixes=ncol(X)
    XX=matrix(1,length(y),1)
    for(var in 1:Fixes) XX=cbind(XX,model.matrix(~X[,var]-1))
    X=XX
    rm(XX,Fixes)
  }
  
  # Defaults of Z: making "NULL","formula" and "matrix" as "list"
  if(is.null(Z)&is.null(iK)) stop("Either Z or iK must be specified")
  if(is.null(Z)) Z=list(diag(length(y)))
  if(class(Z)=="matrix") Z = list(Z)
  if(class(Z)=="formula"){
    Z=model.frame(Z)
    Randoms=ncol(Z)
    ZZ=list()
    for(var in 1:Randoms) ZZ[[var]]=model.matrix(~Z[,var]-1)
    Z=ZZ
    rm(ZZ,Randoms)
  }
  
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
  
  # Some parameters with notation from the book
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
  if(any(is.na(y))){
    MIS = which(is.na(y))
    W=W[-MIS,]
    y=y[-MIS]
    if(!is.null(iR)) iR=iR[-MIS,-MIS]
  }
  n = length(y)
  
  if(is.null(iR)){
    r = crossprod(W,y)
    WW = (crossprod(W))
  }else{
    r = crossprod(W,iR)%*%y
    WW = (crossprod(W,iR))%*%W
  }
  
  # Covariance Matrix
  Sigma = matrix(0,N,N)
  for(i in 1:Randoms) Sigma[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]] = iK[[i]]*lambda[i]
  # Matching WW and Sigma
  C = WW+Sigma
  
  g = rep(0,N)
  
  # Variance components
  dfu = q+DF
  dfe = n+DF
  
  # Saving memory for some vectors
  e = rep(0,N)
  
  # convergence factor
  VC = c(Va,Ve)
  cf = 1
  
  # LOOP
  while(cf>1e-8){
    
    # Ve/Va
    lambda = Ve/Va
    
    # Updating C
    for(i in 1:Randoms) Sigma[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]] = iK[[i]]*lambda[i]
    C = WW+Sigma
    
    # Residual variance
    e = y - tcrossprod(g,W)
    Ve = (tcrossprod(e)+S*DF)/dfe
    
    # the C++ SAMP updates "g" and doesn't return anything
    gs(C,g,r,N)
    
    # Random variance
    for(i in 1:Randoms){
      SS = S*DF+(sum(crossprod(g[Qs1[i+1]:Qs2[i+1]],iK[[i]])*(g[Qs1[i+1]:Qs2[i+1]])))
      SS = max(SS,1e-8)
      Va[i] = SS/dfu[i]
    }
    
    # convergence
    cf = sum((c(Va,Ve)-VC)^2)
    VC = c(Va,Ve)
    
  }
  
  names(g) = paste("b",0:(N-1),sep="")
  for(i in 1:Randoms) names(g)[Qs1[i+1]:Qs2[i+1]] = paste("u",i,".",1:Qs0[i+1],sep="")
  names(VC) = c(paste("Va",1:Randoms,sep=""),"Ve")
  
  # List of Coefficients
  
  RESULTS = list(
    "Coef" = g,
    "VC" = VC,
    "Fit" = W1%*%g
  )
  
  # Return
  return( RESULTS )
  
}

plot.gibbs = function(x,...){
  anyNA = function(x) any(is.na(x))
  par(ask=TRUE)
  vc = nrow(x$Posterior.VC)-1
  for(i in 1:vc) plot(density(x$Posterior.VC[i,],...),main=paste("Posterior: Term",i,"variance"))
  plot(density(x$Posterior.VC[vc+1,]),main=paste("Posterior: Residual Variance"),...)
  par(ask=FALSE)
}

covar = function(sp=NULL,rho=3.5,type=1,dist=2.5){
  if(is.null(sp)) {
    sp = cbind(rep(1,49),rep(c(1:7),7),as.vector(matrix(rep(c(1:7),7),7,7,byrow=T)))
    colnames(sp)=c("block","row","col")
    cat("Example of field information input 'sp'\n")
    print(head(sp,10))
  }
  if(type==1) cat("Exponential Kernel\n")
  if(type==2) cat("Gaussian Kernel\n")
  obs=nrow(sp)
  sp[,2]=sp[,2]*dist
  fields=sp[,1]
  Nfield=length(unique(fields))
  quad=matrix(0,obs,obs)
  for(j in 1:Nfield){
    f=which(fields==j);q=sp[f,2:3]    
    e=dist(q);e=as.matrix(e);e=round(e,3)
    d=exp(-e^type/(rho^type))
    d2=round(d,2);quad[f,f]=d2} 
  if(obs==49){
    M=matrix(quad[,25],7,7)
    dimnames(M) = list(abs(-3:3),abs(-3:3))
    print(M)
  } 
  if(obs!=49) return(quad)
}

PedMat=function(ped=NULL){
  if(is.null(ped)){
    id = 1:11
    dam = c(0,1,1,1,1,2,0,4,6,8,9)
    sire = c(0,0,0,0,0,3,3,5,7,3,10)
    example = cbind(id,dam,sire)
    cat('Example of pedigree\n')
    print(example)
    cat('- It must follow chronological order\n')
    cat('- Zeros are used for unknown\n')
  }else{
    n = nrow(ped)
    A=diag(0,n)
    for(i in 1:n){
      for(j in 1:n){
        if(i>j){A[i,j]=A[j,i]}else{
          d = ped[j,2]
          s = ped[j,3]
          if(d==0){Aid=0}else{Aid=A[i,d]}
          if(s==0){Ais=0}else{Ais=A[i,s]}
          if(d==0|s==0){Asd=0}else{Asd=A[d,s]}
          if(i==j){Aij=1+0.5*Asd}
          if(i!=j){Aij=0.5*(Aid+Ais)}
          A[i,j]=Aij}}}
    return(A)}}
