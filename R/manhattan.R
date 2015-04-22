plot.NAM = function(x,...,alpha=0.05,colA=2,colB=4){
   
  gwas=x
  chr=as.numeric(summary(factor(as.numeric(gwas$MAP[,1]))))
      
  if(gwas$Method=="P3D"){
    FGWASplot=function(Fgwas,chr,AA,BB,alpha,...){
      col=c();for(i in 1:length(chr)){if((i%%2)==0){b=AA}else{b=BB};a=rep(b,chr[i]); col=c(col,a)}
      W = Fgwas$PolyTest$wald
      plot(1:length(W),W,col=col,xlab="Chromosome",ylab="Wald Statistic",xaxt = "n",...)}
    # Plot
    if (is.null(alpha)){FGWASplot(gwas,chr=chr,AA=colA,BB=colB,alpha=alpha,...)}
    else {FGWASplot(gwas,chr=chr,AA=colA,BB=colB,alpha=alpha,...)
          A=1-alpha;LRmax=qchisq(A,0.5);lim=-log(dchisq(LRmax, 0.5));abline(h=lim,col=1,lty=2)}
    
    }else{
    
    RGWASplot=function(Rgwas,chr,AA,BB,alpha,...){
      LRT=(Rgwas$PolyTest$lrt);for(i in 1:length(LRT)){if(LRT[i]<0){LRT[i]=0}}
      if(!is.null(alpha)){
        LRT = -log(dchisq(LRT,df=0.5))
        funLRT = function(lrt){if(lrt<0|lrt==Inf|lrt==-Inf){lrt=0};return(lrt)}
        LRT = sapply(LRT,FUN = funLRT)}
      col=c();for(i in 1:length(chr)){if((i%%2)==0){b=AA}else{b=BB};a=rep(b,chr[i]); col=c(col,a)}
      plot(1:length(Rgwas$PolyTest$lrt),LRT,col=col,xlab="Chromosome",xaxt = "n",...)}
    # Plot
    if (is.null(alpha)){RGWASplot(gwas,chr=chr,AA=colA,BB=colB,alpha=alpha,...)}
    else {RGWASplot(gwas,chr=chr,AA=colA,BB=colB,alpha=alpha,ylab="-log(p-value)",...)
      A=1-alpha;LRmax=qchisq(A,0.5);lim=-log(dchisq(LRmax, 0.5));abline(h=lim,col=1,lty=2)}
    }
  medians=rep(NA,length(chr))
  for(i in 1:length(chr)) medians[i] = median(which(gwas$MAP[,1]==i))
  axis(1, at=round(medians), labels=1:length(medians))
}