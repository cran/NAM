
SSM = function(y,dta,gen,it=300,bi=100,th=1,...){
 
  # dta should include a column called "ID Block Row Col"
  # Pkgs required: Matrix and NAM
  
  cat("Checking for missing and setting up parameters \n")
  # Dealing with missing information
  if(anyNA(y)){
    w=which(y)
    y=y[-w]
    cov=dta[-w,]
  }
  n = length(y)
  gen = data.matrix(gen)
  
  # Missing loci?
  if(anyNA(gen)){
   expectation = function(x){
     if(anyNA(x)) x[is.na(x)] = mean(x)
     return(x)}
   E = apply(gen,2,expectation)
   dim(E) = dim(gen)
   gen = E; rm(E)} 
  
  # Isolating spatial information
  if(sum(c("Block","Row","Col")%in%names(dta))==3){
    spline = TRUE
    spdta = dta[,c("Block","Row","Col")]
  }else{
    spline = FALSE
    cat("No spline: 'dta' does not have columns Block-Row-Col \n")
  }
  
  # Isolating genotypic information
  id = dta[,'ID']
  
  # Isolating of fixed effects
  if(any(!names(dta)%in%c("Block","Row","Col","ID"))){
    fixed = TRUE
    dta2 = dta[,-which(names(dta)%in%c("Block","Row","Col","ID"))]
    dta2 = data.frame(dta2)
    X = Matrix::sparse.model.matrix(~.-1,data=dta2)
  }else{
    fixed = FALSE
    cat("No additional effects other than splines and genomics \n")
  }
  
  # Search for NN plots
  NNsrc = function(sp=NULL,rho=1,dist=3,...){
    # OCTREE search function 
    NN = function(a,sp){
      # shared environment
      s0 = which(sp[,1]==a[1])
      # row neighbors
      s1 = s0[which(abs(sp[s0,2]-a[2])<(rho))]
      s2 = s1[which(abs(sp[s1,3]-a[3])<=(rho*dist))]
      # col neighbors
      s3 = s0[which(abs(sp[s0,3]-a[3])<(rho*dist*0.5))]
      s4 = s3[which(abs(sp[s3,2]-a[2])<=(rho))]
      # Neighbors only
      s5 = union(s2,s4)
      wo = which(sp[s5,1]==a[1]&sp[s5,2]==a[2]&sp[s5,3]==a[3])
      s5 = s5[-wo]
      return(s5)}
    # Set the data format
    rownames(sp) = 1:nrow(sp)
    if(is.data.frame(sp)) sp = data.matrix(sp)
    # Search for neighbors within environment
    Local_NN = apply(sp,1,NN,sp=sp)
    cat(paste('Average n# of components:',round(mean(sapply(Local_NN,length)),2)),'\n')
    # RETURN
    return(Local_NN)
  }
  
  # Sparse matrix of NN
  NNmat = function(NN){
    names(NN) = paste(names(NN),'.',sep='')
    j = unlist(NN)
    i = as.numeric(gsub('\\..+','',names(j)))
    n = tail(i,1)[1]
    x = Matrix::sparseMatrix(i,j,x=1,dims = c(n,n))
    return(x)
  }
  
  # Covariate from NNmat
  SPcov = function(R,e){
    k = Matrix::crossprod(R,e)
    r = Matrix::rowSums(R)
    k = as.vector(k/r)
    k = suppressWarnings(k-mean(k,na.rm=TRUE))
    if(anyNA(k)) k[is.na(k)]=0
    return(k)
  }
  
  if(isTRUE(spline)){
    # Preparing kNN matrix
    cat("Preparing sparse mapping matrix (spline) \n")
    NN = NNsrc(spdta)
    R = NNmat(NN)
  }
  
  # Preparing genotypic mapping matrix
  cat("Generating the incidence matrix of genotypes \n")
  Z = Matrix::sparse.model.matrix(~id-1)
  colnames(Z) = gsub('^id','',colnames(Z))
  Weight = Matrix::colSums(Z)
  
  # Getting the first round fitted
  # (a) Fixed
  fitSM = function(y,X,rdg=1){
    XX = Matrix::crossprod(X)
    Matrix::diag(XX) = Matrix::diag(XX)+rdg
    Xy = Matrix::crossprod(X,y)
    b = Matrix::solve(XX,Xy)
    hat = Matrix::tcrossprod(Matrix::t(b),X)  
    e = y - as.vector(hat)
    ret = list('b'=as.vector(b),'e'=e,'y'=hat)
  }
  # intercept
  mu = mean(y)
  e = y-mu
  # other "fixed" effects with minimal ridge
  if(isTRUE(fixed)){
    LE = fitSM(e,X)
    e = LE$e
    b = LE$b
  }
  
  # (b) Preparing genotypic data
  cat("Computing priors and genomic parameters \n")
  MAP = function(e,Z,Weight) as.vector(Matrix::crossprod(Z,e))/Weight
  xx = apply(gen,2,crossprod)
  p = ncol(gen)
  g = rep(0,p)
  L = rep(1,p)
  Ve = 1
  
  # (c) Running first round
  E = MAP(e,Z,Weight) 
  update = NAM::KMUP(X=gen,b=g,xx=xx,E=E,L=L,p=p,Ve=Ve,pi=0)
  
  # First round of WGR: Setting priors
  R2 = 0.5
  df_prior = 5
  MSx = sum(apply(gen, 2, var, na.rm = T))
  S_prior = R2 * var(y, na.rm = T) * (df_prior + 2)/MSx
  shape_prior = 1.1
  rate_prior = (shape_prior - 1)/S_prior
  
  # BV
  g = update$b
  bv = Matrix::tcrossprod(g,gen)
  bvs = as.vector(Matrix::tcrossprod(Z,bv))
  e = e-bvs
  
  # Update BayesA variance components
  S_conj = rgamma(1, p * df_prior/2 + shape_prior, sum(1/(g^2))/2 + rate_prior)
  # For every round
  Vb = (S_conj + g^2)/rchisq(p, df_prior + 1)
  S_conj = rgamma(1, p * df_prior/2 + shape_prior,sum(1/Vb)/2 + rate_prior)  
  # Update Ve and Lambda
  Ve = crossprod(e)/rchisq(1,n+2)
  L = Ve/Vb
  
  # (d) Adding Splines to X
  if(isTRUE(spline)){
    spEff = SPcov(R,e)
    Q = qr.solve(spEff,e)*spEff
    e = e - Q
    sp = Q
  }
  
  # What to STORE
  MCMC = seq(bi,it,th)
  mc = length(MCMC)
  G = VE = Mu = 0
  VB = rep(0,p)
  if(isTRUE(spline)) SP = rep(0,n) 
  if(isTRUE(fixed)) B = rep(0,ncol(X))
  
  # LOOP STARTS HERE
  cat("Fitting the model \n")
  pb = txtProgressBar(style = 3)
  for(i in 1:it){
    
    # (a) Intercept
    muE = mean(e)
    mu = mu+muE
    e = e-muE
    
    # (b) Other fixed
    if(isTRUE(fixed)){
      LE = fitSM(e,X)
      b = LE$b + b
      e = LE$e
    }
    
    # (c) Update genetic components
    bvs0 = bvs
    E = MAP(e,Z,Weight)
    update = NAM::KMUP(X=gen,b=g,xx=xx,E=E,L=L,p=p,Ve=Ve,pi=0)
    g = update$b
    bv = Matrix::tcrossprod(g,gen)
    bvs = as.vector(Matrix::tcrossprod(Z,bv))
    e = e+bvs0-bvs
    
    # (d) Update VC
    Va = (sum(g^2) + S_prior)/rchisq(1, df_prior + p) # Ridge
    Vb = (S_conj + g^2)/rchisq(p, df_prior + 1) # BayesA
    Vc = abs(g)*mean(abs(g))*0.1 # Bayesian LASSO
    Vm = (Va+Vb+Vc)/3 # Average
    S_conj = rgamma(1, p * df_prior/2 + shape_prior,sum(1/Vb)/2 + rate_prior)  
    Ve = crossprod(e)/rchisq(1,n+2)
    L = Ve/Vm
    
    # (e) Update Splines
    if(isTRUE(spline)){
      spEff = SPcov(R,e)
      Q = qr.solve(spEff,e)*spEff
      e = e - Q
      sp = Q + sp
    }
    
    # (f) Storage
    if(i %in% MCMC){
      Mu = Mu+mu
      G = G+g
      VB = VB+Vb
      VE = VE+Ve
      if(isTRUE(fixed)) B = B+b
      if(isTRUE(spline)) SP = SP+sp
    }
    
    # (g) Update status bar
    setTxtProgressBar(pb, i/it)
  }
  close(pb)
  
  # Posterior means
  Mu = Mu/mc
  G = G/mc
  VB = VB/mc
  VE = VE/mc
  
  # Fitted
  RND = as.vector(Matrix::tcrossprod(Z,Matrix::tcrossprod(G,gen)))
  hat = Mu + RND
  names(RND) = rownames(gen)
  
  if(isTRUE(spline)){
    SP = SP/mc
    hat = Mu + SP
  } 
  
  if(isTRUE(fixed)){
    B = B/mc; names(B)=colnames(X)
    FIX = as.vector(Matrix::tcrossprod(X,B))
    hat = hat + FIX
  }
  
  FINAL = list('hat'=hat,'obs'=y,
               "mu"=Mu,"Z"=Z,"g"=G,
               "Vg"=VB,"Ve"=VE,
               "EBV"=RND)
  
  if(isTRUE(spline)){
    FINAL = c(FINAL,list("sp"=SP))
  }
  
  if(isTRUE(fixed)){
    FINAL = c(FINAL,list("X"=X,"b"=B))
  }
  
  class(FINAL) = "SSM"
  return(FINAL) 
  
}

plot.SSM = function(x,...){
  if("g"%in%names(x)){
    par(mfrow=c(2,1))
    plot((x$g)^2,pch=20,ylab=expression(beta^2),xlab='GENOME',main='Allele Effects')
    plot(x$hat,x$obs,pch=20,xlab="FITTED",ylab="OBSERVED",main='Model Fitness')
    lines(x=c(-1e10,1e10),y=c(-1e10,1e10),type='l',lty=3,lwd=3,col=2)
  }else{
    par(mfrow=c(2,1))
    plot(x$hat,x$obs,pch=20,xlab="FITTED",ylab="OBSERVED",main='Model Fitness')
    lines(x=c(-1e10,1e10),y=c(-1e10,1e10),type='l',lty=3,lwd=3,col=2)
  }
}


