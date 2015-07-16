report = function(y,gen,fam=NULL,chr=NULL,exp=TRUE,reml=FALSE){

  if(reml){
    G = tcrossprod(gen); G = G/mean(diag(G))
    MM = reml( y = y, K = G )
    ZZ = crossprod(gen)
    lambda = as.numeric(MM$VC[2]/MM$VC[1])
    diag(ZZ) = diag(ZZ)+rep(lambda,ncol(gen))
    RHS = crossprod(gen,y-mean(y)) 
    Joint.Effect = solve(ZZ,RHS)
    h2= round(as.numeric(MM$VC[3]),2)
  }else{
    MM = gibbs(y,Z=gen)
    h2 = (MM$VC.estimate[1])/sum(MM$VC.estimate)
    Joint.Effect = as.vector(MM$Coef.estimate[-1])
  }
  
if(anyNA(y)){mis = which(is.na(y));y=y[-mis];
if(!is.null(fam)){fam=fam[-mis]};gen=gen[-mis,]}  

ID = gwas2(y=y,gen=gen,fam=fam,chr=chr)

LRT = ID$PolyTest$lrt
LOD = round(LRT/4.61,2)
P.value = round(-log(dchisq(LRT,0.5)),2)
P.value[P.value<0] = 0

Var.SNP =ID$PolyTest$tau_k
Var.Exp.by.SNP = Var.SNP / var(y,na.rm = T)

Vx = apply(gen,2,var,na.rm=T)
COV = cov(y,gen)
R2 = t(COV^2)
Marg.Effect = t(COV/Vx)

r1 = data.frame(ID$SNPs,ID$MAP[,1:2],LOD,
        P.value,Var.Exp.by.SNP,R2,Marg.Effect,Joint.Effect)

if(!is.null(fam)){u = unique(fam);f = length(u)
  for(i in 1:f){
    w = which(fam==u[i])
    eff = t(cov(y[w],gen[w,])/Vx)
    r1 = cbind(r1,eff)
names(r1)[8+i] = paste("Eff.Pop.",u[i],sep='')}}

if(exp){
  write.csv(r1,"Report.GWAS.csv",row.names=F)
  leg = paste("(h2 = ",h2,")",sep='')
  png("Report.GWAS.png",width = 800, height = 500)
  plot(ID,main="Genome-wide Association Studies")
  legend("topleft",leg,bty = "n")
  dev.off()}

LIST=list(GWAS=ID,MixMod=MM,Report=r1)
return(LIST)
}