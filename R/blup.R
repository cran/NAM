
# Average Information REML

blup = function(y,Z=NULL,X=NULL,K=NULL,Z2=NULL,Z3=NULL,K2=NULL,K3=NULL,conv=1e-8,maxit=50,verb=TRUE){

  if(is.null(K)&is.null(Z))stop("Either Z or K must be specified")
    
  Z1=Z
  K1=K

  tr = function(x){tr=sum(diag(x));return(tr)}
  n = length(y)
  I = matrix(0, ncol=2, nrow=2);  s = matrix(0, ncol=1, nrow=2)
  
  
# Allowing formulas for fixed variables
  if(class(X)=="formula"){X=model.matrix(X)}else{
    if(class(X)=="matrix"){X=X}else{
    X=matrix(1,n,1)}}

# Test X for Singularity
if(any(eigen(crossprod(X))$values<=0)) stop("Fixed effects ran into a singularity.")

# Allowing formulas for random variables
  if(is.null(Z1)){}else{if(class(Z1)=="formula"){Z1=model.matrix(Z1);Z1[,1]=Z1[,1]*(rowSums(Z1)==1)}}
  if(is.null(Z2)){}else{if(class(Z2)=="formula"){Z2=model.matrix(Z2);Z2[,1]=Z2[,1]*(rowSums(Z2)==1)}}
  if(is.null(Z3)){}else{if(class(Z3)=="formula"){Z3=model.matrix(Z3);Z3[,1]=Z3[,1]*(rowSums(Z3)==1)}}

# cell.means = model.matrix(~as.factor(Z)-1)

# Checking if dimensions are compatible
  if(nrow(X)!=length(y))stop("Response variable and fixed effect matrix have incompatible dimensions")
  if(is.null(Z1)){}else{if(nrow(Z1)!=n){stop("Response variable and Z1 have incompatible dimensions")}}
  if(is.null(Z2)){}else{if(nrow(Z2)!=n){stop("Response variable and Z2 have incompatible dimensions")}}
  if(is.null(Z3)){}else{if(nrow(Z3)!=n){stop("Response variable and Z3 have incompatible dimensions")}}

# Checking the existence of kinships
  if(is.null(Z1)){}else{if(is.null(K1)){K1=diag(ncol(Z1))}}
  if(is.null(Z2)){}else{if(is.null(K2)){K2=diag(ncol(Z2))}}
  if(is.null(Z3)){}else{if(is.null(K3)){K3=diag(ncol(Z3))}}

# Checking the existence of design matrices
  if(is.null(K1)){}else{if(is.null(Z1)){Z1=diag(n)}}
  if(is.null(K2)){}else{if(is.null(Z2)){Z2=diag(n)}}
  if(is.null(K3)){}else{if(is.null(Z3)){Z3=diag(n)}}

# Accomodating missing values
# This step must be done after generating Kinships

if(anyNA(y)){
  
  missing = which(is.na(y))
  y.old = y
  y = y[-missing]
  n2 = length(y)
  X.old = X
  X = as.matrix(X[-missing,])
  
  if(is.null(Z1)){}else{Z1.old = Z1}
  if(is.null(Z2)){}else{Z2.old = Z2}
  if(is.null(Z3)){}else{Z3.old = Z3}
  
  if(is.null(Z1)){}else{Z1 = Z1[-missing,]}
  if(is.null(Z2)){}else{Z2 = Z2[-missing,]}
  if(is.null(Z3)){}else{Z3 = Z3[-missing,]}
}else{
  
  y.old = y
  n2 = length(y)
  X.old = X
  if(is.null(Z1)){}else{Z1.old = Z1}
  if(is.null(Z2)){}else{Z2.old = Z2}
  if(is.null(Z3)){}else{Z3.old = Z3}
    
}


# Final number of random variables
  RanVar = 4-(is.null(Z1)+is.null(Z2)+is.null(Z3))

  # Eigendecomposition step
if(RanVar==2){cat("Eigendecomposing Kinship",'\n') # Eigen decomposing K1
  eig.A = eigen(K1,symmetric=TRUE);  D1 = eig.A$values; D1[D1<=0]=10^-8; U1 = eig.A$vectors}

# Starting AIREML
cat("Iterating Average Information REML",'\n')

# Numerators
V1=Z1%*%(K1)%*%t(Z1); 
if(is.null(Z2)){}else{V2=Z2%*%(K2)%*%t(Z2)}
if(is.null(Z3)){}else{V3=Z3%*%(K3)%*%t(Z3)} 

# AI matrix and s updater
I = matrix(0, ncol=RanVar, nrow=RanVar); 
s = matrix(0, ncol=1, nrow=RanVar)

nVar = rep(var(y)*3,RanVar)
nVar[RanVar]=nVar[RanVar]/2
Var = matrix(nVar,ncol=1)
diffVar = rep(1,RanVar)

MaxDif = 1
iter = 0
while(MaxDif>conv && iter<maxit){  

iter = iter+1
  
    lambda1 = Var[1]/Var[RanVar]
    if(is.null(Z2)){}else{lambda2 = Var[2]/Var[RanVar]}
    if(is.null(Z3)){}else{lambda3 = Var[3]/Var[RanVar]}
    
    # Solution of inverse V
    if(RanVar==2){Vinv = (Z1 %*% U1%*%chol2inv(diag(D1*lambda1+1))%*%t(U1)%*%t(Z1))/Var[RanVar]
                  }else{
      V = diag(n2)*Var[RanVar] + Z1%*%(K1*Var[1])%*%t(Z1)
      if(is.null(Z2)){}else{V = V + Z2%*%(K2*Var[2])%*%t(Z2)}
      if(is.null(Z3)){}else{V = V + Z3%*%(K3*Var[3])%*%t(Z3)}
      Vinv=chol2inv(V)
    }
    
    # Updating the P matrix    
    if(ncol(X)==1){P=Vinv}else(P=Vinv-Vinv%*%X%*%chol2inv(t(X)%*%Vinv%*%X)%*%t(X)%*%Vinv)
    
    # Updating the Average Information matrix
    PV1 = P%*%V1; 
    if(is.null(Z2)){}else{PV2 = P%*%V2}
    if(is.null(Z3)){}else{PV3 = P%*%V3} 
    
    yP=t(y)%*%P;  Py=P%*%y;
    
    if(RanVar==4){
    # Average-Information Matrix for 3 random
    I[1,1] = tr( yP%*%PV1%*%PV1%*%Py ); I[1,2] = tr( yP%*%PV1%*%PV2%*%Py ); I[1,3] = tr( yP%*%PV1%*%PV3%*%Py ); I[1,4] = tr( yP%*%PV1%*%Py );
    I[2,1] = I[1,2]; I[2,2] = tr( yP%*%PV2%*%PV2%*%Py ); I[1,3] = tr( yP%*%PV2%*%PV3%*%Py ); I[2,4] = tr( yP%*%PV2%*%Py );
    I[3,1] = I[1,3]; I[3,2] = I[2,3]; I[3,3] = tr( yP%*%PV3%*%PV3%*%Py ); I[2,4] = tr( yP%*%PV3%*%Py );
    I[4,1] = I[1,4]; I[4,2] = I[2,4]; I[4,3] = I[3,4]; I[4,4] = tr(yP%*%Py)
    }else{if(RanVar==3){
    # Average-Information Matrix for 2 random
    I[1,1] = tr( yP%*%PV1%*%PV1%*%Py ); I[1,2] = tr( yP%*%PV1%*%PV2%*%Py ); I[1,3] = tr( yP%*%PV1%*%Py );
    I[2,1] = I[1,2]; I[2,2] = tr( yP%*%PV2%*%PV2%*%Py ); I[2,3] = tr( yP%*%PV2%*%Py );
    I[3,1] = I[1,3]; I[3,2] = I[2,3]; I[3,3] = tr(yP%*%Py)
    }else{
    # Average-Information Matrix for 1 random
    I[1,1] = tr( yP%*%PV1%*%PV1%*%Py ); I[1,2] = tr( yP%*%PV1%*%Py );
    I[2,1] = I[1,2]; I[2,2] = tr(yP%*%Py)}}
    
    AI=0.5*I;
    diag(AI)[diag(AI)<=0]=10e-8
    
    s[1,1] = -0.5* (tr(PV1)-(t(y)%*%PV1%*%Py ))
    if(is.null(Z2)){}else{s[2,1] = -0.5*(tr(PV2)-(t(y)%*%PV2%*%Py ))}
    if(is.null(Z3)){}else{s[3,1] = -0.5*(tr(PV3)-(t(y)%*%PV2%*%Py ))}
    s[RanVar,1] = -0.5*(tr(P)-(t(y)%*%P%*%Py ))
    
    # UPDATING VARIANCE:  v* = AI^-1 x s
    newVar = Var + chol2inv(AI)%*%s
    newVar[newVar<=0]=conv*100
    newVar[is.nan(newVar)]=conv*100
    newVar[is.infinite(newVar)]=conv*100
    
    # Checking convergence
    diffVar = abs(Var - newVar)
    MaxDif = max(diffVar)
    
    Var = newVar

if(isTRUE(verb)){
  cat("Iter =",iter,"Conv =",max(diffVar),
      paste(c(paste("ICC",1:(RanVar-1)," = ",sep="")),
            round(Var[1:(RanVar-1)]/sum(Var),3)),'\n')}

}

VarComp = Var
rownames(VarComp)=paste(c(paste("Var",1:(RanVar-1),sep=""),"VarE"))

# SOLVE MIXED MODEL
cat("Solving Henderson Model",'\n')

if(RanVar==4){
    
  XX=crossprod(X);      Z1X=crossprod(Z1,X);   Z2X=crossprod(Z2,X);  Z3X=crossprod(Z3,X);
  XZ1=crossprod(X,Z1);  Z1Z1=crossprod(Z1)+chol2inv(K1)*lambda1;     Z2Z1=crossprod(Z2,Z1); Z3Z1=crossprod(Z3,Z1);
  XZ2=crossprod(X,Z2);  Z1Z2=crossprod(Z1,Z2); Z2Z2=crossprod(Z2)+chol2inv(K2)*lambda2;   Z3Z2=crossprod(Z3,Z2);
  XZ3=crossprod(X,Z3);  Z1Z3=crossprod(Z1,Z3); Z2Z3=crossprod(Z2,Z3); Z3Z3=crossprod(Z3,Z3)+chol2inv(K3)*lambda3;
  
  col1=rbind(XX,Z1X,Z2X,Z3X)
  col2=rbind(XZ1,Z1Z1,Z2Z1,Z3Z1)
  col3=rbind(XZ2,Z1Z2,Z2Z2,Z3Z2)
  col4=rbind(XZ3,Z1Z3,Z2Z3,Z3Z3)
  
  LHS = cbind(col1,col2,col3,col4)
  RHS = rbind( t(X)%*%y, t(Z1)%*%y, t(Z2)%*%y, t(Z3)%*%y )
  LHS = LHS + diag(1e-4,ncol(LHS))
  W = qr.solve(LHS,RHS,tol=1e-10)
  
  b = W[1:ncol(X)]; W=W[-c(1:ncol(X))]
  u1 = W[1:ncol(Z1)]; W=W[-c(1:ncol(Z1))]
  u2 = W[1:ncol(Z2)]; W=W[-c(1:ncol(Z2))]
  u3 = W
  Fit = X.old%*%b + Z1.old%*%u1 + Z2.old%*%u2 + Z3.old%*%u3
  BLUP = list("VarComp"=VarComp[,1],"b"=b,"u1"=u1,
              "u2"=u2,"u3"=u3,"Fit"=Fit,"Obs"=y.old)
  
  
}else{if(RanVar==3){
  
  XX=crossprod(X);      Z1X=crossprod(Z1,X);   Z2X=crossprod(Z2,X);
  XZ1=crossprod(X,Z1);  Z1Z1=crossprod(Z1)+chol2inv(K1)*lambda1;     Z2Z1=crossprod(Z2,Z1); 
  XZ2=crossprod(X,Z2);  Z1Z2=crossprod(Z1,Z2); Z2Z2=crossprod(Z2)+chol2inv(K2)*lambda2;
  
  row1=cbind(XX,XZ1,XZ2)
  row2=cbind(Z1X,Z1Z1,Z1Z2)
  row3=cbind(Z2X,Z2Z1,Z2Z2)
  
  LHS = rbind(row1,row2,row3)
  LHS = LHS + diag(1e-4,ncol(LHS))
  RHS = rbind( t(X)%*%y, t(Z1)%*%y, t(Z2)%*%y )
  W = qr.solve(LHS,RHS,tol=1e-10)
  
  b = W[1:ncol(X)]; W=W[-c(1:ncol(X))]
  u1 = W[1:ncol(Z1)]; W=W[-c(1:ncol(Z1))]
  u2 = W
  Fit = X.old%*%b + Z1.old%*%u1 + Z2.old%*%u2

  BLUP = list("VarComp"=VarComp[,1],"b"=b,"u1"=u1, "u2"=u2,"Fit"=Fit,"Obs"=y.old)
  
}else{

XX=crossprod(X);      Z1X=crossprod(Z1,X);
XZ1=crossprod(X,Z1);  Z1Z1=crossprod(Z1)+chol2inv(K1)*lambda1;

row1=cbind(XX,XZ1)
row2=cbind(Z1X,Z1Z1)

LHS = rbind(row1,row2)
LHS = LHS + diag(1e-4,ncol(LHS))
RHS = rbind( t(X)%*%y, t(Z1)%*%y )
W = qr.solve(LHS,RHS,tol=1e-10)

b = W[1:ncol(X)]
u=W[-c(1:ncol(X))]
Fit = X.old%*%b + Z1.old%*%u

BLUP = list("VarComp"=VarComp[,1],"b"=b,"u"=u,"Fit"=Fit,"Obs"=y.old)
}}

return(BLUP)
}