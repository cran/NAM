snpQC=function(gen,psy=1,MAF=0.05,remove=TRUE,impute=FALSE){
  anyNA = function(x) any(is.na(x))
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
rf <- function(xmis){ # Author: D.Stekhoven, stekhoven@stat.math.ethz.ch
  maxiter = 10; ntree = 100
  mtry = floor(sqrt(ncol(xmis)))
  cutoff=NULL;classwt=NULL;strata=NULL;replace=T
  sampsize=NULL;nodesize=NULL;maxnodes=NULL
  n <- nrow(xmis);  p <- ncol(xmis)
  ## perform initial S.W.A.G. on xmis (mean imputation)
  ximp <- xmis; xAttrib <- lapply(xmis, attributes)
  varType <- character(p);  for (t.co in 1:p){
    if (is.null(xAttrib[[t.co]])){ varType[t.co] <- 'numeric'
      ximp[is.na(xmis[,t.co]),t.co] <- mean(xmis[,t.co], na.rm = TRUE)
    }else{ varType[t.co] <- 'factor'
     max.level <- max(table(ximp[,t.co]))
     class.assign <- sample(names(which(max.level == summary(ximp[,t.co]))), 1)
     if (class.assign != "NA's"){ximp[is.na(xmis[,t.co]),t.co] <- class.assign } else {
     while(class.assign=="NA's"){class.assign=sample(names(which(max.level==summary(ximp[,t.co]))), 1)}
     ximp[is.na(xmis[,t.co]),t.co] <- class.assign }}}
  NAloc <- is.na(xmis);  noNAvar <- apply(NAloc, 2, sum) 
  sort.j <- order(noNAvar); sort.noNAvar <- noNAvar[sort.j]
  nzsort.j <- sort.j[sort.noNAvar > 0]; Ximp <- vector('list', maxiter)
  iter <- 0; k <- length(unique(varType));  convNew <- rep(0, k)
  convOld <- rep(Inf, k); OOBerror <- numeric(p); names(OOBerror) <- varType
  if (k==1){if(unique(varType)=='numeric'){names(convNew)=c('numeric')}else{names(convNew)=c('factor')}
    convergence <- c(); OOBerr <- numeric(1)} else {
    names(convNew) <- c('numeric', 'factor');convergence <- matrix(NA, ncol = 2);OOBerr <- numeric(2)}
   stopCriterion <- function(varType, convNew, convOld, iter, maxiter){
    k <- length(unique(varType)); if (k == 1){(convNew < convOld) & (iter < maxiter)} else {
    ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & (iter < maxiter)}}
  while (stopCriterion(varType, convNew, convOld, iter, maxiter)){
    if (iter != 0){convOld <- convNew;OOBerrOld <- OOBerr}
    cat("RF iteration", iter+1, "\n");t.start <- proc.time();ximp.old <- ximp
    for(s in 1:p){varInd=sort.j[s];if(noNAvar[[varInd]]!=0){obsi=!NAloc[,varInd];misi=NAloc[, varInd];
     obsY=ximp[obsi,varInd];obsX=ximp[obsi, seq(1,p)[-varInd]];misX=ximp[misi, seq(1,p)[-varInd]];
     typeY=varType[varInd];RF=randomForest(x=obsX,y=obsY,ntree=ntree,mtry=mtry,replace=T,
     sampsize=if(!is.null(sampsize))sampsize[[varInd]] else if (replace) nrow(obsX) else
       ceiling(0.632*nrow(obsX)), nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
       maxnodes = if (!is.null(maxnodes)) maxnodes else NULL); OOBerror[varInd] <- RF$mse[ntree]
       misY <- predict(RF, misX); ximp[misi, varInd] <- misY }}; iter <- iter+1;Ximp[[iter]] <- ximp
    t.co2 <- 1;for (t.type in names(convNew)){ t.ind <- which(varType == t.type)
    convNew[t.co2]=sum((ximp[,t.ind]-ximp.old[,t.ind])^2)/sum(ximp[,t.ind]^2);t.co2=t.co2 + 1}}
  if (iter == maxiter){out <- Ximp[[iter]]}else{out <- Ximp[[iter-1]]}
  return(out)}
     gen=gen3;gen[gen==5]=NA   
     k=100*length(which(is.na(gen)))/length(gen)
     k=round(k,2);cat(k,"% of missing data",'\n')
     cat("Imputations being performed by Random Forest",'\n')
    if(any(is.na(gen))){gen=suppressWarnings(rf(gen));gen=round(gen)}
    gen3=gen}
return(gen3)}

# function to remove repeated genotypes
cleanREP = function(y,fam,gen,thr=0.95){
  if(is.vector(y)) y=matrix(y,ncol=1)
  GG=function(gen,r=1){
    a1=(gen-1)
    a1[a1==-1]=0
    A1=(tcrossprod(a1))
    a2=-(gen-1)
    a2[a2==-1]=0
    A2=(tcrossprod(a2))
    d=round(exp(-abs(gen-1)))
    D=tcrossprod(d)
    G=A1+A2+D;G=(G/ncol(gen))^r
    return(G)}
  cat('solving identity matrix\n')
  G = GG(gen) # identity
  rownames(G) = 1:nrow(G)
  lt = G*lower.tri(G) # lower triang
  r = 1* lt>thr # logical matrix: repeatitions
  # starting point of new data
  rownames(gen) = 1:nrow(gen)
  Ny=y;  Nfam=fam;  Ngen=gen
  # summary
  cs = colSums(r) # how many times id will be repeated
  while(any(cs>0)){
    i = which(cs>0)[1]
    cat("indiviual",rownames(gen)[i],"had",cs[i],'duplicate(s)\n')
    w = which(r[,i])
    if(ncol(y)>1){y[i,] = colMeans(y[c(i,w),],na.rm=T)
    }else{y[i] = mean(y[c(i,w)],na.rm=T)}
    if(ncol(y)>1){Ny=Ny[-w,]}else{Ny=Ny[-w]}
    Nfam=Nfam[-w]
    Ngen=Ngen[-w,]
    r = r[-w,]
    cs = colSums(r)
  }
  return(list(y=Ny,gen=Ngen,fam=Nfam))
}
