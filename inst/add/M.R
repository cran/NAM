
M = function(...,RUN=FALSE){
  
  # Load packages
  require(Rcpp)
  require(Matrix)
  require(compiler)
  
  # Load C++
  cppFunction(depends = 'RcppEigen',code =
  'SEXP GS2Eigen(Eigen::Map<Eigen::VectorXd> e, Eigen::MappedSparseMatrix<double> X,
              Eigen::Map<Eigen::VectorXd> b, Eigen::MappedSparseMatrix<double> XX,double Lmb) {
  int P = X.cols();
  int N = X.rows();
  Eigen::VectorXd Y(N);
  Eigen::VectorXd r(P);
  Y = X * b + e;
  r = X.transpose() * Y;
  double b0;
  Eigen::VectorXd Xi;
  for(int i=0; i<P; i++){
  b0 = b(i);
  Xi = XX.col(i);
  b(i) = ( r(i) - Xi.transpose()*b + Xi(i)*b0  ) / (Xi(i)+Lmb);
  }
  e = Y - X * b;
  return List::create(Named("b")=b,Named("e")=e);
  }')
  
  cppFunction(depends = 'RcppEigen',code ='
  SEXP GSEigen(Eigen::Map<Eigen::VectorXd> e,
  Eigen::MappedSparseMatrix<double> X,
  Eigen::Map<Eigen::VectorXd> b,
  Eigen::Map<Eigen::VectorXd> XX,
  double Lambda) {
  int N = X.cols();
  double oldB;
  Eigen::VectorXd Xi;
  for(int i=0; i<N; i++){
  oldB = b(i);
  Xi = X.col(i);
  b(i) = ( Xi.transpose() * e + XX(i) * oldB ) / (XX(i)+Lambda);
  e = e - Xi * (b(i)-oldB);
  }
  return List::create(Named("b")=b,Named("e")=e);
  }')
  
  cppFunction(code ='
  SEXP GSANA(NumericVector e, NumericMatrix gen, NumericVector b, NumericVector Lmb,
  NumericVector xx, double cxx, double Ve, double h2 = 0.5, int maxit = 75){
  double tol = 10e-10;
  double phi = cxx*((1-h2)/h2);
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  double b0,eM;
  double mu = mean(e);
  e = e-mu;
  // Regulation coefficients
  NumericVector Vb(p);
  double b1;
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
  // Regression coefficients loop
  bc = b+0;
  for(int j=0; j<p; j++){
  // Ordinary Least Square
  b0 = b[j];
  b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(Lmb(j)+xx(j));
  b[j] = b1;
  // Residuals update
  e = e-gen(_,j)*(b1-b0);}
  // Intercept update
  eM = mean(e);
  mu = mu+eM;
  e = e-eM;
  // Variance components
  Vb = b*b+(Ve/(xx+Lmb));
  Lmb = sqrt(phi*Ve/Vb);
  // Convergence
  ++numit;
  cnv = sum(abs(bc-b));
  if( cnv<tol ){break;}
  }
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b); }
  // Output
  return List::create(Named("mu")=mu,
  Named("b")=b,
  Named("e")=e,
  Named("h")=fit,
  Named("Lmb")=Lmb,
  Named("h2")=sum(Vb)/(sum(Vb)+Ve));
  }')
  
  cppFunction(code = '
  SEXP GSRUK(NumericVector e, NumericMatrix gen, NumericVector b,
  NumericVector eV, double cxx, double h2 = 0.5){
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  double eM;
  double mu = mean(e);
  e = e-mu;
  NumericVector Lmb = 1 + cxx*((1-h2)/h2)/(eV);
  double G, b0;
  for(int j=0; j<p; j++){
  // Ordinary Least Square
  b0 = b[j];
  G = (sum(gen(_,j)*e)+b0)/(Lmb[j]); 
  b[j] = G;
  // Residuals update
  e = e-gen(_,j)*(b[j]-b0);
  eM = mean(e);
  mu = mu+eM;
  e = e-eM;
  }
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b); }
  // Output
  return List::create(Named("mu")=mu,Named("b")=b,Named("e")=e,Named("h")=fit);
  }')
  
  cppFunction(code = '
  SEXP MSX(NumericMatrix X){
  int p = X.ncol(); int n = X.nrow(); double m; NumericVector xx(p); NumericVector sx(p);
  for(int k=0; k<p; k++){ xx[k] = sum(X(_,k)*X(_,k)); m = sum(X(_,k)); sx[k] = m*m/n; }
  double cxx = sum(xx-sx)/(n-1); return List::create(Named("cxx")=cxx,Named("xx")=xx);}')
  
  cppFunction(code = '
  void IMP(NumericMatrix X){;int p = X.ncol(); int n = X.nrow();
  LogicalVector MIS(n); NumericVector x(n); NumericVector z; double EXP;
  for(int j=0; j<p; j++){;if(is_true(any(is_na(X(_,j))))){
  x = X(_,j); MIS = is_na(x);z = x[!MIS]; EXP = mean(z);
  X(_,j) = ifelse(MIS,EXP,x);};};}')
  
  # Load and compile main function
  eMM = function(y,random=NULL,fixed=NULL,data=NULL,gen=NULL,ped=NULL,binomial=FALSE,QC=FALSE,R2=0.5,maxit=NULL,verb=TRUE){
    
    # Base algorithm settings
    if(!is.null(data)) data = droplevels.data.frame(data); as = 0
    it = ifelse(is.null(maxit),yes =  ifelse(is.null(gen),20,10),no = maxit)
    if(!verb){ cat = function(a) a }
    
    ###########################
    ### Data Quality Control ##
    ###########################
    
    if(QC&!is.null(data)){
      
      cat("Performing quality control\n")
      
      # QC on fixed (factor levels <3)
      if(!is.null(fixed)){
        a = attr(terms(fixed),"term.labels")
        w = c()
        # 3 iterations
        for(j in 1:3 ){
          for(i in a){
            if(is.factor(data[[i]])){
              b = table(data[[i]]); if(any(b<3)){
                w = names(which(b<3)); w = which(data[[i]] %in% w)
                data = droplevels.data.frame(data[-w,])}
            }else if(is.numeric(data[[i]])){
              if(any(is.na(data[[i]]))){
                w = which(is.na(data[[i]]));
                data = droplevels.data.frame(data[-w,])}}}}
        rm(i,j,a,b,w)}
      
      # QC on random (factor levels <3)
      if(!is.null(random)){
        a = attr(terms(random),"term.labels")
        w = c()
        # 3 iterations
        for(j in 1:3){
          for(i in a){
            if(is.factor(data[[i]])){
              b = table(data[[i]]); if(any(b<3)){
                w = names(which(b<3)); w = which(data[[i]] %in% w)
                data = droplevels.data.frame(data[-w,])}
            }else if(is.numeric(data[[i]])){
              if(any(is.na(data[[i]]))){
                w = which(is.na(data[[i]]));
                data = droplevels.data.frame(data[-w,])}}}}
        rm(i,j,a,b,w)}
      
      # QC on y (outliers)
      Y = data[[deparse(substitute(y))]]
      sdy = sqrt(var(Y,na.rm = T))
      my = mean(Y,na.rm = T)
      upLim = my+3*sdy
      loLim = my-3*sdy
      w = which(Y<loLim | Y>upLim)
      if(length(w)>0) data = droplevels.data.frame(data[-w,])
      rm(Y,sdy,my,upLim,loLim,w)
    }
    
    ######################
    ### Design Matrices ##
    ######################
    
    # Response variable
    if(!is.null(data)) y = data[[deparse(substitute(y))]]
    
    # Fixed effects or intercept
    if(!is.null(fixed)){
      cat("Setting fixed effect design matrices\n")
      if(anyNA(data[, attr(terms(fixed),"term.labels") ])) stop("Missing parameters are not allowed")
      fixed = update(fixed,~.-1)
      X = sparse.model.matrix(fixed,data = data)
      XX = as(crossprod(X),"dgCMatrix")
      B = rep(0,ncol(X)); names(B) = colnames(X)
      Lmb_X = (1e-12)*sum(colMeans(X^2)-colMeans(X)^2)
      mu = 0
      FIX = TRUE
    }else{
      mu = 0
      FIX = FALSE
    }
    
    # Random effects
    if(!is.null(random)){
      cat("Setting random effect design matrices\n")
      rnd = attr(terms(random),"term.labels")
      nr = length(rnd)
      Z = list()
      if(anyNA(data[,rnd])) stop("Missing parameters are not allowed")
      for(i in 1:nr) Z[[i]] = sparse.model.matrix(formula(paste('~',rnd[i],'-1')),data = data)
      names(Z) = rnd
      U = lapply(Z, function(x){z=rep(0,ncol(x));names(z)=colnames(x);return(z)}  )
      RND = TRUE
    }else{
      nr = 0
      rnd = NULL
      RND = FALSE
    }
    
    # Copy of Z of ID for correct prediction later
    if( RND & 'ID'%in%rnd ){
      Z_ID = Z$ID
      Zblup = TRUE
    }else{
      Zblup = FALSE
    }
    
    # Missing values
    if(anyNA(y)){
      y0 = y;
      if(FIX) X0 = X
      if(RND) Z0 = Z
      wNA = which(is.na(y))
      y = y[-wNA];
      if(FIX) X = X[-wNA,]
      if(RND) Z = lapply(Z,function(z) z[-wNA,] )
      MIS = TRUE
    }else{
      MIS = FALSE
    }
    n = length(y)
    
    # Binomial
    if(binomial){
      cat("Logit tranformation\n")
      ry = range(y)
      y = (y-min(y))/(max(y-min(y)))
      y = round(y,3)
      y = ifelse(y==0,0.001,y)
      y = ifelse(y==1,0.999,y)
      y = log(y/(1-y))
    }
    
    # Random parameters for regularization
    if(RND){
      ZZ = lapply(Z, function(z) as(crossprod(z),"dgCMatrix") ) 
      Z_cxx = lapply(Z, function(X) sum(colMeans(X^2)-colMeans(X)^2) )
      Lmb_Z = mapply( function(cxx,h2) cxx*((1-h2)/h2), cxx = Z_cxx, h2=R2)
      if('ID'%in%rnd){ZXX = diag(ZZ$ID); names(ZXX)=sub('ID','',colnames(ZZ$ID)) } 
      trAC22 = sapply(ZZ,function(x) mean(1/(diag(x)+1)))
      fxx = ifelse(FIX,sum(1/diag(XX)),0)
      WCW = (sum(unlist(mapply(function(ZZ,lmb)sum(diag(ZZ)/(diag(ZZ)+lmb)),ZZ=ZZ,lmb=as.list(Lmb_Z))))+fxx)/n
    }
    
    # Genotypes
    if(!is.null(gen)){
      cat("Computing genomic parameters")
      
      # Controls
      if( RND & is.null(rownames(gen)) ) stop("Matrix 'gen' does not have row names to indentify genotypes")
      if( RND & !"ID"%in%rnd ) stop("No random effect called 'ID' was declared")
      if(RND){
        colnames(Z$ID) = gsub("^ID","",colnames(Z$ID))
        names(U$ID) = gsub("^ID","",names(U$ID))
        proportion_genotyped = round(mean(colnames(Z$ID)%in%rownames(gen))*100,2)
        cat(' (',proportion_genotyped," percent of IDs are genotyped)\n",sep='')
        if(proportion_genotyped==0) stop("IDs from data and genotypic matrix do not match\n")
      }
      
      # Imputations
      if(!is.matrix(gen)) gen = data.matrix(gen)
      if(anyNA(gen)) IMP(gen)
      
      # Mapping matrix
      if(RND){
        wGen = which(!rownames(gen)%in%colnames(Z$ID))
        if(length(wGen)!=0) gen = gen[-wGen,]
        gFit = rep(0,ncol(Z$ID)); names(gFit) = colnames(Z$ID)
        gHat = rep(0,n)
      }
      
      # Parameters
      msx = MSX(gen)
      M_xx = msx$xx
      M_cxx = msx$cxx
      Beta = rep(0,ncol(gen))
      LMB = rep(M_cxx,ncol(gen))
      
      # Logical for Genomic
      GEN = TRUE }else{ GEN = FALSE }
    
    # Index individuals without genotype
    if( GEN & 'ID'%in%rnd ){
      gebv = U$ID
      wg = which( colnames(Z$ID) %in% rownames(gen) )
      ng = names(gebv)[ which(names(gebv) %in% rownames(gen)) ]
      ZID = Z$ID
      Z$ID[,wg] = 0
      ALL_GENOTYPED = GEN & (ncol(Z$ID)==length(wg)) 
    }else{ ALL_GENOTYPED = FALSE; wg=0 }
    
    # SibRidge function that returns an incidence matrix N, so that NN' = A
    SibRidge = function(ped){
      # formatting missing
      ped[is.na(ped)] = 0
      no = nrow(ped)
      # formatting into character
      if(is.data.frame(ped)) ped = apply(ped,2,as.character)
      # get levels
      k = data.matrix(ped); k = sort(unique(as.vector(k)))
      k = k[k!='0']
      np = length(k)
      # set sparse
      x = sqrt( 1 - 0.25 * apply(ped[,2:3],1,function(x) sum(x!=0)) )
      SP = data.frame(I=1:no,J=NA,X=x)
      SP$J = as.numeric(factor(ped[,1],levels = k))
      N = sparseMatrix(i=SP$I,j=SP$J,x=SP$X,dims = c(no,np),dimnames = list(ped[,1],k))
      # add parentals
      for(i in 2:3){
        SP$J = as.numeric(factor(ped[,i],levels = k))
        SP$X = 0.5
        w = which(!is.na(SP$J))
        if(length(w)>0) N = N + sparseMatrix(i=SP$I[w],j=SP$J[w],x=SP$X[w],dims = c(no,np),dimnames = list(ped[,1],k))
      }
      # Finale
      return(N)
    }
    
    # Pedigree
    if( (!is.null(ped)) & (!ALL_GENOTYPED) ){
      cat("Computing pedigree parameters")
      
      pen = SibRidge(ped)
      
      # Controls
      if( RND & !"ID"%in%rnd ) stop("No random effect called 'ID' was declared")
      if(RND){
        colnames(Z$ID) = gsub("^ID","",colnames(Z$ID))
        names(U$ID) = gsub("^ID","",names(U$ID))
        proportion_ped = round(mean(colnames(Z$ID)%in%rownames(pen))*100,2)
        cat(' (',proportion_ped," percent of IDs have pedigree)\n",sep='')
        if(proportion_ped==0) stop("IDs from data and pedigree do not match\n")
      }
      
      # Mapping matrix
      if(RND){
        wPen = which(!rownames(pen)%in%colnames(Z$ID))
        if(length(wPen)!=0) pen = pen[-wGen,]
        pFit = rep(0,ncol(Z$ID)); names(pFit) = colnames(Z$ID)
        pHat = rep(0,n)
      }
      
      # Parameters
      P_xx = apply(pen,2,function(X) sum(X^2))
      P_cxx = sum(apply(pen,2,var))
      pBeta = rep(0,ncol(pen))
      
      PED = TRUE
      
    }else{ PED = FALSE }
    
    # Index individuals without pedigree
    if( PED & 'ID'%in%rnd ){
      pebv = U$ID
      wp = which( colnames(Z$ID) %in% rownames(pen) )
      np = names(pebv)[ which(names(pebv) %in% rownames(pen)) ]
      pZID = Z$ID
      Z$ID[,wp] = 0
    }
    
    # Starting value for variances
    Se = Ve = var(y)*(1-R2)
    if(RND){ Va = rep(var(y)*R2,nr); names(Va) = rnd } 
    e = y
    
    # Remove overall mean
    m0 = mean(e)
    mu = mu+m0
    e = e-m0
    
    ###########
    ### Loop ##
    ###########
    
    if(verb) pb = txtProgressBar(style = 3)
    for(i in 1:it){
      
      # Fixed coefficient update
      if(FIX) upFix = GS2Eigen(e,X,B,XX,Lmb_X)
      
      # Randomc oefficient update
      if(RND){
        for(j in 1:nr){
          # Special case for ID (in case of purely WGR)
          if( GEN & rnd[j] =='ID' ){
            # If all individuals are genotyped, it skips.
            if(ncol(Z$ID)!=length(wg)){
              upRnd = GS2Eigen(e,Z[[j]],U[[j]],ZZ[[j]],Lmb_Z[j])
            }
          }else{
            upRnd = GS2Eigen(e,Z[[j]],U[[j]],ZZ[[j]],Lmb_Z[j])
          }
        }}
      
      # Genomic coefficient update
      if(GEN){
        if(RND){
          bvs0 = gHat
          Gmap = rowSums(crossprod(ZID,e)/ZXX)[rownames(gen)]
          G_update = GSANA(Gmap,gen,Beta,LMB,M_xx,M_cxx,Ve,R2,40)
          ggFit = ( Gmap*ZXX[rownames(gen)] + G_update$h*Lmb_Z["ID"] ) / ( ZXX[rownames(gen)]+ Lmb_Z["ID"] )
          gFit = U$ID
          gFit[rownames(gen)] = ggFit
          U$ID = gFit
          # Done
          gHat[1:n] = (ZID%*%gFit)[,1]
          e = e + bvs0 - gHat
        }else{
          G_update = GSANA(e,gen,Beta,LMB,M_xx,M_cxx,Ve,R2,40)
        }
      }
      
      # Pedigree coefficient update
      if( PED & (!ALL_GENOTYPED) ){
        if(RND){
          bvs0 = pHat
          Pmap = rowSums(crossprod(pZID,e)/ZXX)[rownames(pen)]
          P_update = GSEigen(Pmap,pen,pBeta,P_xx,P_cxx*(1-R2)/R2);
          pFit[rownames(pen)] = as.vector(pen%*%pBeta)
          pHat[1:n] = (pZID%*%pFit)[,1] 
          e = e + bvs0 - pHat
        }else{P_update = GSEigen(e,pen,pBeta,P_xx,P_cxx*(1-R2)/R2)}}
      
      # Update overall intercept
      m0 = mean(e)
      mu = mu+m0
      e = e-m0
      
      # Variance components update
      if(RND){
        #  Merging GV with EBV and GEBV (where GEBV override EBV)
        if( GEN & 'ID'%in%rnd ){ U$ID = gFit  }
        if( PED & 'ID'%in%rnd & (!ALL_GENOTYPED) ){
          U$ID[rownames(pen)] = P_update$h }
        
        # Update variance components
        trAC22 = (unlist(mapply(function(ZZ,lmb)sum(1/(diag(ZZ)+lmb)),ZZ=ZZ,lmb=as.list(Lmb_Z))))/sapply(U,length)
        Ve = crossprod(e,y)[1,1]/(n+ifelse(FIX,ncol(X),1))
        Va = sapply(U,var) + trAC22*Ve
        Lmb_Z = (Ve/Va)
        
      }
      
      # Update progress bar
      if(verb) setTxtProgressBar(pb, i/it)
    }
    
    if(verb) close(pb)
    
    ##################
    ### Wrapping up ##
    ##################
    
    # Add fixed effects
    if(FIX){
      fx = list( Fixed = list(X=X,B=B,B0=mu))
    }else{
      fx = list( Intercept = mu )
    }
    
    # Add random effects
    if(Zblup) Z$ID = Z_ID
    if(RND){
      fx[[2]] = list(Z=Z,U=U)
      fx[[3]] = list(Va=round(Va,3),Ve=round(Ve,3))
      names(fx)[2:3] = c('Random','VarComp')
    }
    
    if(GEN){
      names(Beta) = colnames(gen)
      fx[[ifelse(RND,4,2)]] = Beta
      names(fx)[ifelse(RND,4,2)] = paste("MarkerEffects",round(G_update$h2,2)*100,sep='.')
    }
    
    if(PED){
      fx[[length(fx)+1]] = list(Ped_Z=pen,Ped_u=pBeta)
      names(fx)[length(fx)] = "Ped"
    }
    
    # Fitting model
    if(MIS){
      fit = rep(0,length(y0))
      if(FIX){fit=fit+X0%*%B+mu}else{fit=fit+mu}
      if(RND) for(i in 1:nr) fit=fit+Z0[[i]]%*%U[[i]]
      fx[[(length(fx)+1)]] = as.vector(fit)
      names(fx)[length(fx)] = "Hat"
      if(binomial){
        b_fit = exp(fit)/(1+exp(fit));
        b_fit = b_fit*(ry[2]-ry[1])+ry[1]
        fx[[(length(fx)+1)]] = round(as.numeric(b_fit))
        names(fx)[length(fx)] = "Binomial_Hat"
      } 
    }else{
      fit = rep(0,n)
      if(FIX){fit=fit+X%*%B+mu}else{fit=fit+mu}
      if(RND) for(i in 1:nr) fit=fit+Z[[i]]%*%U[[i]]
      if(GEN&!RND) fit=fit+gen%*%Beta
      fx[[(length(fx)+1)]] = as.vector(fit)
      names(fx)[length(fx)] = "Hat"
      if(binomial){
        b_fit = exp(fit)/(1+exp(fit));
        b_fit = b_fit*(ry[2]-ry[1])+ry[1]
        fx[[(length(fx)+1)]] = round(as.numeric(b_fit))
        names(fx)[length(fx)] = "Binomial_Hat"
      }
    }
    
    # Log-likelihood, R2, Adj R2, BIC
    loglik = sum(dnorm(as.numeric(y),as.numeric(fit),as.numeric(Ve),T))/2
    px = 1+ifelse(FIX,ncol(X),0)+ifelse(RND,length(Z),0)+ifelse(GEN|PED,1,0)
    if(MIS){
      R2 = cor(as.numeric(y0),as.numeric(fit),use='p')^2
    }else{
      R2 = cor(as.numeric(y),as.numeric(fit),use='p')^2
    }
    adj_R2 = 1 - ( (1-R2)*(n-1) / (n-px-1) )
    BIC = log(n)*px + log(Ve)*n
    if(MIS){ obs = y0 }else{ obs = y } 
    stat = list(loglik=loglik,R2=R2,adj_R2=adj_R2,BIC=BIC,obs=obs)
    # AdjR2
    fx[[(length(fx)+1)]] = stat
    names(fx)[length(fx)] = "stat"
    
    # After debuggin, recompile this pack with the following line of code
    # require(roxygen2);require(devtools);roxygenise('eMM');use_rcpp('eMM');document('eMM');build('eMM');check('eMM')
    # install.packages("eMM_1.2.5.tar.gz", repos = NULL, type = "source")
    
    class(fx) = "eMM"
    return(fx)
  }
  cmpfun(eMM)
  
  # End
  if(RUN){
    return(eMM(...))
  }else{
    cat("Arguments (RUN = TRUE): \nM(y, random = NULL, fixed = NULL, data = NULL,\n   gen = NULL, ped = NULL, binomial = FALSE,\n   QC = FALSE, R2 = 0.5, maxit = NULL, verb = TRUE)\n")
  } 
  
}
