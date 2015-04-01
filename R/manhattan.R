manhattan <-
  function(gwas,colA=2,colB=4,alpha=0.05,...){
    
    chr=as.numeric(summary(factor(as.numeric(gwas$MAP[,1]))))
    
    RGWASplot=function(Rgwas,chr,AA,BB,alpha,...){
      
      LRT=(Rgwas$PolyTest$lrt);for(i in 1:length(LRT)){if(LRT[i]<0){LRT[i]=0}}
      
      if(is.null(alpha)==FALSE){
        LRT = -log(dchisq(LRT,df=0.5))
        #for(i in 1:length(LRT)){if(LRT[i]==Inf|LRT[i]==-Inf|LRT[i]<0){LRT[i]=0}}
        funLRT = function(lrt){if(lrt<0|lrt==Inf|lrt==-Inf){lrt=0};return(lrt)}
        LRT = sapply(LRT,FUN = funLRT)
      }
      
      col=c();for(i in 1:length(chr)){if((i%%2)==0){b=AA}else{b=BB};a=rep(b,chr[i]); col=c(col,a)}
      
      plot(1:length(Rgwas$PolyTest$lrt),LRT,col=col,xlab="Genome",...)}
    
    
    if(is.null(alpha)){
      RGWASplot(gwas,chr=chr,AA=colA,BB=colB,alpha=alpha, ...)
    }else{
      RGWASplot(gwas,chr=chr,AA=colA,BB=colB,alpha=alpha,ylab="-log(p-value)", ...)
      A = 1-alpha
      LRmax = qchisq(A,0.5)
      lim = -log(dchisq(LRmax,0.5))
      abline(h=lim,col=1,lty=2)
    }
  
  }