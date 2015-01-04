Fst=function(gen,fam){
    m=ncol(gen)
    INDEX=fam
    FUN=function(allele){ #Fst
      AA=length(which(allele==2))
      Aa=length(which(allele==1))
      aa=length(which(allele==0))
      PA=(AA+0.5*Aa)/(AA+Aa+aa)
      Hn=2*PA*(1-PA)
      return(c(PA,Hn))}
    Fst=c()
    for(i in 1:m){
      pop=tapply(gen[,i],INDEX,FUN)
      pop2=matrix(unlist(pop),ncol=2,byrow=T)
      Hs=mean(pop2[,2])
      Ht=2*mean(pop2[,1])*(1-mean(pop2[,1]))
      st=1-Hs/Ht
      if(is.nan(st)) st=0
      if(is.na(st)) st=0
      if(is.infinite(st)) st=0
      Fst=c(Fst,st)}
    
    C = 1-0.1*(4+6*abs(Fst)/max(abs(Fst)))
    
    plot(1:length(Fst),Fst,
      main="Wright's Fst",pch=20,xlab="Genome",
      col=(rgb(C,C,C)))
    
    return(Fst)}