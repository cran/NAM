reml<-function(y,X=NULL,Z=NULL,K=NULL){

  N = length(y) 
# Dealing with fixed effect matrix
  if(is.null(X)){X = matrix(1,N,1)}else{
    if(is.matrix(X)){if(nrow(X)!=N) stop("Fixed effect does not match dimensions of response variable")
    }else{if(class(X)=="formula"){X=model.matrix(X)}}}

# Dealing with random effect
if(is.null(K)&is.null(Z))stop("Either Z or K must be specified")

if(is.null(K)){
  if(class(Z)=="formula"){Z=model.matrix(Z)-1}
  V=tcrossprod(Z)}

if(is.null(Z)){V=K}

if(is.null(Z)!=T&&is.null(K)!=T){
  if(class(Z)=="formula"){Z=model.matrix(Z)-1}
  V=crossprod(Z,K);V=tcrossprod(V,Z)}

K=V

# Function starts here
  m = which(is.na(y)) # missing values
  if(anyNA(y)){
    y=y[-m];x=X[-m,]
    k=K[m,-m];K=K[-m,-m]
  }else{x=X}

# Defining log-REML
  loglike<-function(theta){
    lambda<-exp(theta)
    logdt<-sum(log(lambda*delta+1))
    h<-1/(lambda*delta+1)
    yy<-sum(yu*h*yu)
    yx<-matrix(0,q,1)
    xx<-matrix(0,q,q)
    for(i in 1:q){
      yx[i]<-sum(yu*h*xu[,i])
      for(j in 1:q){xx[i,j]=sum(xu[,i]*h*xu[,j])}}
    loglike = -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
    return(-loglike)}
  fixed<-function(lambda){
    h<-1/(lambda*delta+1)
    yy<-sum(yu*h*yu)
    yx=timesVec(yu,h,xu,q)
    xx=timesMatrix(xu,h,xu,q,q)
    beta<-qr.solve(xx,yx)
    sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
    sigma2 = as.numeric(sigma2)
    var<-diag((chol2inv(xx))*sigma2)
    stderr<-sqrt(var)
    return(c(beta,stderr,sigma2))}

# Eigendecomposition of K
  qq<-eigen(as.matrix(K),symmetric=T)
  delta<-qq[[1]]
  uu<-qq[[2]]
  q<-max(ncol(x),1)
  n<-ncol(K)
  vp<-var(y)
  yu<-t(uu)%*%y
  xu<-t(uu)%*%x
  theta<-0
 
# Finding lambda through optimization
   parm<-optim(par=theta,fn=loglike,method="L-BFGS-B",lower=-10,upper=10)
   lambda<-exp(parm$par)
 
# More efficient optimazer
# require(mcGlobaloptim);lambda=exp(multiStartoptim(theta,loglike,method="nlminb")$par)

# Variance components and fixed effect coefficient
  parmfix<-fixed(lambda)
  beta<-parmfix[1:q]
  sd<-parmfix[(q+1):(2*q)]
  B = cbind(beta,sd)
  Ve<-parmfix[2*q+1]; Vg<-lambda*Ve
  h2=Vg/(Vg+Ve); VC = data.frame(Vg,Ve,h2)

# Random effect coefficient and prediction
# VanRaden(2008)
  re =(y-x%*%beta)
  iG = timesMatrix(uu,1/(lambda*delta+1),uu,n,n)
  U = iG%*%re
  if(length(m)>0){
    C=k%*%iG%*%U
    hat = rep(0,N)
    hat[m] = C
    hat[-m] = U
    U = hat}
return(list("VC"=VC,"Fixed"=B,"EBV"=U))}