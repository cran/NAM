manhattan <-
  function(gwas,colA=2,colB=4,GenDist=FALSE,OtherDist=NULL,...){
    
    chr=as.numeric(summary(factor(as.numeric(gwas$MAP[,1]))))
    
    RGWASplot=function(Rgwas,chr,AA,BB,GD,OD,...){
      
      LRT=(Rgwas$PolyTest$lrt);for(i in 1:length(LRT)){if(LRT[i]<0){LRT[i]=0}}
      
      col=c();for(i in 1:length(chr)){if((i%%2)==0){b=AA}else{b=BB};a=rep(b,chr[i]); col=c(col,a)}
      
      if(is.null(OD)!=TRUE){plot(OD,LRT,col=col,xlab="Genome",...)}
      
      if(GD==TRUE){plot(Rgwas$MAP[,6],LRT,col=col,xlab="Genome",...)}
      
      plot(1:length(Rgwas$PolyTest$lrt),LRT,col=col,xlab="Genome",...)}
    
    RGWASplot(gwas,chr=chr,AA=colA,BB=colB,GD=GenDist,OD=OtherDist, ...)
  }