snpQC=function(gen,psy=1,MAF=0.05,remove=TRUE,impute=FALSE){
  # CHECKING REDUNDANT MARKERS
  gen2=gen; redundancy=c(0); for(i in 1:(ncol(gen)-1)){
    a=mean((gen[,i]==gen[,(i+1)]),na.rm=TRUE);redundancy=c(redundancy,a)} 
  a=which(redundancy>psy);b=length(which(redundancy>psy))
  if(b>0){cat("Genotypic data contains",b,"redundant SNPs",'\n')
          if(remove==TRUE){gen2=gen[,-a]}
  }else{cat("No redundant SNPs found",'\n')}
  # CHECKING MINOR ALLELE FREQUENCY
  if(MAF>0){
    LAF=c();for(j in 1:(ncol(gen2))){
      AA=length(which(gen2[,j]==2))
      Aa=length(which(gen2[,j]==1))
      aa=length(which(gen2[,j]==0))
      Total=AA+Aa+aa
      PA=(AA+0.5*Aa)/Total
      Pa=(aa+0.5*Aa)/Total
      lowerAF=min(PA,Pa)
      LAF=c(LAF,lowerAF)}
    maf=which(LAF<MAF)
    hist(LAF,col=3,nclass=50,main="Histogram of MAF",xlab="Minor Allele Frequency")
    if(length(maf)>0){
          cat("There are",length(maf),"markers with MAF below the threshold",'\n')
          if(remove==TRUE){gen3=gen2[,-maf]}
      }else{cat("No marker below MAF threshold",'\n');gen3=gen2}
  }else{gen3=gen2}
  if(impute){
     gen=gen3;gen[gen==5]=NA   
     k=100*length(which(is.na(gen)))/length(gen)
     k=round(k,2);cat(k,"% of missing data",'\n')
     cat("Imputations being performed by Random Forest",'\n')
    # IMPUTATION by Random Forest implemented in the missForest package
    if(anyNA(gen)){options(warn=-1);gen=missForest(gen);gen=gen$ximp
        gen=round(gen,0);options(warn=0)}
    gen3=gen}
return(gen3)}