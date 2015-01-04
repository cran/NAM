reference <-
function(gen,ref=NULL){
  if(is.null(ref)){
    maf=function(z){
      Z=z
      z[z==2]="A"; nA=length(which(z=="A"))
      z[z==0]="B"; nB=length(which(z=="B"))
      if(nA==nB){x=Z}
      if(nA>nB){z[z=="A"]=2;z[z=="B"]=0;x=z}
      if(nA<nB){z[z=="B"]=2;z[z=="A"]=0;x=z}
      return(x)}
    W=apply(gen,2,maf);W=(as.numeric(W));W=matrix(W,ncol=ncol(gen))
    colnames(W)=colnames(gen);
  }else{
    if(ncol(gen)!=length(ref)) stop("Reference parent and matrix of genotypes display non compatible dimensions")
    if(anyNA(ref)||length(which(ref==5))>0) {"Reference parent must have no missing values"}
    gen[is.na(gen)]=5
    change=function(x,y=ref){
      z=c();
      for(i in 1:length(y)){
        if(x[i]==y[i]){I=2}
        else if(x[i]==1){I=1}
        else if(x[i]==5){I=5}
        else if(x[i]!=y[i]){I=0}
        else stop("Matrix must be coded as 012, using 5 or NA for missing data")
        z=c(z,I)}
      return(z)}
    W=apply(gen,1,change);W=t(W);colnames(W)=colnames(gen);
  }
    
  return(W)}
