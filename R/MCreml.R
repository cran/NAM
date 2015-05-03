
MCreml = function(y,K,X=NULL,MC=300,samp=300){
  anyNA = function(x) any(is.na(x))
  if(samp>=length(y)){stop("Sample size has to be smaller than sample space")}
  if(nrow(K)!=length(y)){stop("Kinship and response variable have incompatible dimensions")}
  if(ncol(K)!=length(y)){stop("Kinship and response variable have incompatible dimensions")}
  n = MC; t = samp
  moda=function (x){
    it=5;ny=length(x);k=ceiling(ny/2)-1; while(it>1){
      y=sort(x); inf=y[1:(ny-k)]; sup=y[(k+1):ny]
      diffs=sup-inf; i=min(which(diffs==min(diffs)))
      M=median(y[i:(i+k)]); it=it-1}; return(M)}
  h2 = c(); Vg = c(); Ve = c()
  for(i in 1:n){
    R = sample(1:length(y),t)
    if(any(is.na(y[R]))){mis = which(is.na(y[R])); R=R[-mis] }
    fit = reml(y[R],X=X,K=K[R,R])
    Vg = c( fit$VC[1], Vg )
    Ve = c( fit$VC[2], Ve )
    h2 = c( fit$VC[3], h2 )
  }
  H = unlist(h2); G = unlist(Vg); E = unlist(Ve)
  samples = cbind(G,E,H); rownames(samples) = 1:MC
  mode.Vg = moda(G); mode.Ve = moda(E); mode.h2 = moda(H)
  VC = c(mode.Vg,mode.Ve,mode.h2); names(VC) = c("Vg","Ve","h2")
  result = list("samples"=samples,"modes"=VC)
  
  return(result)
  
}