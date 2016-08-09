# Function to run GenABEL with NAM inputs
GRAMMAR_gamma=function(y,gen,fam=NULL,chr=NULL,cov=NULL){
  if(is.null(chr)) chr = ncol(gen)
  if(any(is.na(gen))) gen = markov(gen,chr)
  # Input GenABEL
  id=paste("A",1:length(y),sep="")
  sex=rep(1,length(y))
  trait=y
  pheno=cbind(id,sex,trait)
  write.table(pheno,"pheno",col.names=T,quote=F,row.names=F)
  name=colnames(gen)
  CHR=c();for(i in 1:length(chr)){a=rep(i,chr[i]);CHR=c(CHR,a)}
  pos=(1:ncol(gen))*100+1
  GEN=t(gen);GEN[GEN==2]="GG";GEN[GEN==1]="GC";
  GEN[GEN==0]="CC";GEN[GEN==5]="00";GEN[is.na(GEN)]="00";
  colnames(GEN)=paste("A",1:length(y),sep="")
  geno=cbind(name,CHR,pos,GEN)
  colnames(geno)[2]="chr"
  write.table(geno,"data.illu",col.names=T,quote=F,row.names=F,sep="\t")
  # GenABLE step
  require(GenABEL)
  convert.snp.illumina(inf="data.illu",out="data.raw",strand="u")
  srdta <- load.gwaa.data(phe="pheno",gen="data.raw",force=TRUE)
  # Kinship
  gm = function(gen,fam) {
    g1 = tcrossprod(gen)
    g1 = g1/mean(diag(g1))
    g2 = ((gen-1)*-1)+1
    g2 = tcrossprod(g2)
    g2 = g2/mean(diag(g2))
    for (i in unique(fam))
      {nam = which(fam == i)
      g1[nam,nam]=g2[nam,nam]+g1[nam,nam]}
    lambda = mean(diag(g1))
    G = g1/mean(diag(g1))
    return(G)}
  if(is.null(fam)){
    G=gen%*%t(gen);G=G/mean(diag(G))
  }else G = gm(gen,fam)
  # GWAS step
  if(is.null(cov)){trait=y}else{trait=y~cov}
  POLY=polygenic(trait,data=srdta,kinship.matrix=G)
  LM=grammar(POLY,data=srdta,method="gamma")
  plot(LM,main="GWAS GRAMMAR-gamma",type="h",lwd=3)
  final = list('LM'=LM,'srdta'=srdta,'POLY'=POLY)
  return(srdta)
}