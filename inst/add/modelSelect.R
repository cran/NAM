#' @title Genomic prediction model evaluation
#' @name Cross-Validation
#' @aliases CV.model CV.check
#' @description Provided the uncertainty associated to the models selection for genomic prediction purposes, this function evaluate a set of model with various properties (mostly with additive nature). With the current setups, this function works reasonably well and sufficiently fast for datasets up to 4-8K SNPs and 3-5k Observations.
#' @usage CV.model(y,gen,k=5,Seeds=1:5,cores=NULL)
#' CV.check(cv)
#' @param y Numeric vector of phenotypes of length  \eqn{n} (missing values are allowed)
#' @param gen Numeric matrix with dimension \eqn{n} x \eqn{p} (missing values are NOT allowed)
#' @param k Indicates the k-fold cross validation, that defines in how many pieces the data will be splitted (5 by default)
#' @param Seeds Seeds being set to pick prediction and validation set. By default is uses 5 seeds: 1,2,3,4,5. The length of 'Seeds' defines how many times the CV will be performed.
#' @param cores Number of cores for parallel processing. If NULL (default), runs in one core.
#' @export CV.model CV.check
#' @return The \code{CV.test} function returns a list of cross-validations (ie. training set and prediction set don't share observations). The additional function \code{CV.check} uses the output of \code{CV.model} to return a list prediction parameters (from best to worst) consisting of:
#' \item{PA}{Prediction accuracy, Pearson correlation between predited and observed values}
#' \item{Rank}{Rank order correlation, Spearman correlation between predited and observed values}
#' \item{MSPE}{Mean square prediction error, standard parameter for supervised machine learning experiments}
#' \item{Bias}{Prediction bias as the slope between predicted and observed values, whre 1 means Unbiased}
#' \item{Top20}{Accuracy (Percentage) of correctly selecting the top 20\% phenotypes}
#' \item{Model_choice}{The five models that had the best average performance for this prediction scenario}
#' @author Alencar Xavier (\url{http://alenxav.wix.com/home})
#' @details The models under evaluation include:
#' 
#' 01) \code{Random Forest}{ (pkg. randomForest) with 500 trees containing as many as 1000 parameters}
#' 
#' 02) \code{Radial Basis Function}{ (pkg. kernlab), with tunning parameter found through cross-validation}
#' 
#' 03) \code{BayesB}{ (pkg. BGLR), Marker effects follow a mixture of binomial/t distribution, with a probability of SNP included around 0.01 sampled from beta-binomial, allows large effect QTLs and null effect markers}
#'
#' 04) \code{BayesC}{ (pkg. BGLR), Marker effects follow a mixture of binomial/normal distribution, with a probability of SNP included around 0.01, sampled from beta-binomial, vide Kuo and Mallick 1998}
#'
#' 05) \code{GBLUP/REML}{ (pkg. rrBLUP), Kernel version of GBLUP, using: ZZ'/alpha, where Z is gen with centrilized allele conding, and alpha is the mean diagonal of ZZ'}
#'
#' 06) \code{RKHS with 1 Gaussian kernel}{ (pkg. BGLR), KERNEL: exp(-E2/median(E2)) where E2 = as.matrix(dist(gen))^2}
#'
#' 07) \code{RKHS with 3 Gaussina kernel}{ (pkg. BGLR), KERNELS: exp(-E2/quantile(E2,0.05)*5), exp(-E2/quantile(E2,0.05)), exp(-E2/quantile(E2,0.05)*0.2)}
#'
#' 08) \code{BayesB + RKHS}{ (pkg. BGLR) Ensemble learning strategy to capture non-linear relationships and large effect QTLs}
#'
#' 09) \code{BayesB + RKHS (3G)}{ (pkg. BGLR) - Ensemble learning strategy, same as above but utilizing three Gaussian kernels}
#'
#' 10) \code{Ridge Regression}{ (pkg. glmnet), where lambda parameter is found through cross-validation}
#'
#' 11) \code{LASSO}{ (pkg. glmnet), where lambda parameter is found through cross-validation}
#'
#' 12) \code{Elastic Net}{ (pkg. glmnet), where lambda parameter is found through cross-validation and alpha parameter set as 0.5}
#'
#' 13) \code{Bayesian Ridge Regression}{ (pkg. BGLR), which simply is a random effect model}
#'
#' 14) \code{Bayesian LASSO}{ (pkg. BGLR) from Park and Casella 2008 (Laplacian shrinkage on L2 loss)}
#'
#' 15) \code{Bayesian Elastic-Net}{ (pkg. bWGR) using bagging and acceptance-rejection to find alpha}
#'
#' 16) \code{Bagging Bayesian Ridge Regression}{ (pkg. bWGR), with 30\% of data used at a time: as a mixture of Ridge and Random forest}
#'
#' 17) \code{Bagging Bayesian Ridge Regression}{ (pkg. bWGR), with 70\% of data used at a time: as a mixture of Ridge and Random forest}
#'
#' 18) \code{k-Nearest Neighbors} { function with k=20, ie., average of 20 most similar individuals based on Euclidean distance}
#'
#' 19) \code{Boosting Regression} { (pkg. gbm) Linear regression with adaptive Boosting}
#'
#' 20) \code{Partial least square regression} { (pkg. pls) with cross-validation to define the number of Eigenvectos}
#'
#' 21) \code{General Ensemble Learning Algorithm} { (pkg. bWGR) GELA: Mix elements of random forest, RKHS, ridge regression and kNN}
#'
#' 22) \code{Weighted General Ensemble Learning Algorithm}{ (pkg. bWGR) wGELA: same as above, but with weighted models based on out-of-bag prediction error}
#'
#' 23) \code{Variational BayesB}{ (pkg. VIGoR) with prior tunning of hyperparameters, starting with Mvar=0.5 and Kappa=0.1 }
#'
#' 24) \code{Variational BayesC}{ (pkg. VIGoR) with the same tunning as in varBayesB }
#'
#' 25) \code{Variational Bayesian stochastic search variable selection}{ (pkg. VIGoR) with the same tunning as in varBayesB }
#'
#' 26) \code{Variational Extended Bayesian LASSO}{ (pkg. VIGoR) with the same tunning as in varBayesB }
#'
#' 27) \code{Variational Bayesian mixture of models}{ (pkg. VIGoR) with the same tunning as in varBayesB }
#'
#' 28) \code{EM-BayesA}{ (pkg. bWGR) with default settings }
#'
#' 29) \code{EM-BayesB}{ (pkg. bWGR) with default settings }
#'
#' 30) \code{EM-BayesC}{ (pkg. bWGR) with default settings }
#'
#' 31) \code{EM-BRR}{ (pkg. bWGR) with default settings }
#'

CV.model=function(y,gen,k=5,Seeds=1:5,cores=NULL){
  
  if(anyNA(gen)) stop("Missing values in 'gen' are not allowed")
  
  CV.test=function(y,gen,k,Seeds){
    # Check if you have packages
    check = function(mypkg) if(!is.element(mypkg,installed.packages()[,1])) install.packages(mypkg)
    check('BGLR')
    check('glmnet')
    check('gbm')
    check('rrBLUP')
    check('pls')
    check('randomForest')
    check('kernlab')
    check('bWGR')
    check('VIGoR')
    # Loading packages
    require(BGLR,quietly = T)
    require(randomForest,quietly = T)
    require(glmnet,quietly = T)
    require(gbm,quietly = T)
    require(rrBLUP,quietly = T)
    require(pls,quietly = T)
    require(kernlab,quietly = T)
    require(bWGR,quietly = T)
    require(VIGoR,quietly = T)
    # KNN function
    knn = function(y,E,k=20){d=function(d) mean(r[which(d<=(sort(d)[k]))]);
    w=which(is.na(y));D=E[,-w];r=y[-w];return(apply(D,1,d))}
    # Kernels
    E2 = as.matrix(dist(gen)^2)
    K=exp(-(E2/mean(E2)))
    K1=exp(-(E2/(quantile(E2,0.05))))
    K2=exp(-(E2/(0.2*quantile(E2,0.05))))
    K3=exp(-(E2/(5*quantile(E2,0.05))))
    eK = eigen(K,symmetric=T)
    eK1 = eigen(K1,symmetric=T)
    eK2 = eigen(K2,symmetric=T)
    eK3 = eigen(K3,symmetric=T)
    Z = apply(gen,2,function(x) x-mean(x))
    G = tcrossprod(Z)
    G = G/mean(diag(G))
    Y = y
    N = nrow(gen)
    cvs = length(Seeds) 
    cat('DONE with EIGENDECOMPOSITIONS\n')
    # Cross-validation function
    folds = function(Seed){
      set.seed(Seed)
      Nk = round(N/k)
      w=sample(1:N,Nk)
      y[w]=NA
      # RF
      f1=randomForest(x=gen[-w,],y=y[-w], mtry = min( floor(ncol(gen)/3), 1000 ))
      cat('RF\n')
      # RBF
      cat('now RBF\n')
      f2=predict(rvm(gen,y),gen)[,1]
      # Bayesian
      f2b = BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BayesB')),verbose=F)
      f2c = BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BayesC')),verbose=F)
      # RKHS
      f3 = BGLR(rmExistingFiles = F,y,ETA=list(list(V=eK$vectors,d=eK$values,model='RKHS')),verbose=F)
      f3b = BGLR(rmExistingFiles = F,y,ETA=list(list(V=eK1$vectors,d=eK1$values,model='RKHS'),
                                                list(V=eK2$vectors,d=eK2$values,model='RKHS'),
                                                list(V=eK3$vectors,d=eK3$values,model='RKHS')),verbose=F)
      cat('Alphabet and HS\n')
      # Lasso
      cv1 = cv.glmnet(x=gen[-w,],y=y[-w],alpha=1)
      lmb1 = cv1$lambda.min
      f4 = glmnet(x=gen[-w,],y=y[-w],lambda = lmb1,alpha = 1)
      f4b = BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BL')),verbose=F)
      # Ridge
      cv2 = cv.glmnet(x=gen[-w,],y=y[-w],alpha=1)
      lmb2 = cv2$lambda.min
      f5 = glmnet(x=gen[-w,],y=y[-w],lambda = lmb2,alpha = 0)
      f5b = BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BRR')),verbose=F)
      # BB+K
      f6 = BGLR(rmExistingFiles = F,y,ETA=list(list(V=eK$vectors,d=eK$values,model='RKHS'),
                                               list(X=gen,model='BayesB')),verbose=F)
      f6b = BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BayesB'),
                                                list(V=eK1$vectors,d=eK1$values,model='RKHS'),
                                                list(V=eK2$vectors,d=eK2$values,model='RKHS'),
                                                list(V=eK3$vectors,d=eK3$values,model='RKHS')),verbose=F)
      cat('done with L1 and L2\n')
      # kNN
      f7 = knn(y,E2,20)
      # Boosting
      f8=gbm::gbm.fit(x=gen[-w,],y=y[-w],distribution="gaussian",n.trees = 150)
      # Bag Bayes
      f9b=wgr(y,gen,it=1800,bi=0,bag=.3,verb=T,rp=FALSE)
      f9c=wgr(y,gen,it=1800,bi=0,bag=.7,verb=T,rp=FALSE)
      # REML GBLUP
      f10 = mixed.solve(y,K=G)
      cat('done with Bags and GBLUP\n')
      # GELA
      f11 = GELA(y,x=gen,weighted = FALSE)
      f12 = GELA(y,gen,weighted = TRUE)
      cat('done with GELAs\n')
      # BagMAN
      f13b = BEN(y=y,gen=gen)
      # EN
      cv3 = cv.glmnet(x=gen[-w,],y=y[-w],alpha=0.5)
      lmb3 = cv3$lambda.min
      f14 = glmnet(x=gen[-w,],y=y[-w],lambda = lmb1,alpha = 0.5)
      # PLS
      cat('now PLS\n')
      f16a=plsr(y[-w]~gen[-w,],validation='CV')
      ncomp=which.max(cor(f16a$validation$pred[,1,],y[-w]))
      pls16a = predict(f16a,gen[w,])[,1,ncomp]
      # varBayes
      f17a = gen[w,] %*% vigor(y,gen,'BayesB',hyperpara(gen,.5,'BayesB',.1),Function = 'tuning')$Beta
      f17b = gen[w,] %*% vigor(y,gen,'BayesC',hyperpara(gen,.5,'BayesC',.1),Function = 'tuning')$Beta
      f17c = gen[w,] %*% vigor(y,gen,'SSVS',hyperpara(gen,.5,'SSVS',.1),Function = 'tuning')$Beta
      f17d = gen[w,] %*% vigor(y,gen,'EBL',hyperpara(gen,.5,'EBL',.1),Function = 'tuning')$Beta
      f17e = gen[w,] %*% vigor(y,gen,'MIX',hyperpara(gen,.5,'MIX',.1),Function = 'tuning')$Beta
      # emBayes
      f18a = gen[w,] %*% emBA(y[-w],gen[-w,])$b
      f18b = gen[w,] %*% emBB(y[-w],gen[-w,])$b
      f18c = gen[w,] %*% emBC(y[-w],gen[-w,])$b
      f18d = gen[w,] %*% emRR(y[-w],gen[-w,])$b
      cat('Done with emBayes and varBayes\n')
      
      # Removing junk files and creating MODs
      file.remove(list.files(pattern='dat'))
      #mod = rep(NA,1+length(grep('^f.+',ls())))
      NamesMod = c('RandomForest','RBF','BayesB','BayesC','RKHS',
                   'RKHS(3GK)','LASSO','BL','RidgeReg',
                   'BRR','BayesB+RKHS','BayesB+3GK','kNN',
                   'Boost','emBRR','Bagging(30%)','Bagging(70%)',
                   'GBLUP-REML','GELA','wGELA','BEN','EN','PLS',
                   'varBayesB','varBayesC','SSVS','EBL','varMIX',
                   'emBayesA','emBayesB','emBayesC','OBSERVATION')
      M = matrix(NA,Nk,length(NamesMod))
      colnames(M) = NamesMod
      # VALUES
      M[,1]=predict(f1,gen[w,]) # RF
      M[,2]=f2[w] # RBF
      M[,3]=f2b$yHat[w] # BayesB
      M[,4]=f2c$yHat[w] # BayesC
      M[,5]=f3$yHat[w] # RKHS
      M[,6]=f3b$yHat[w] # RKHS - 3GK
      M[,7]=predict(f4,gen[w,]) # LASSO
      M[,8]=f4b$yHat[w] # BL
      M[,9]=predict(f5,gen[w,]) # Ridge
      M[,10]=f5b$yHat[w] # BRR
      M[,11]=f6$yHat[w] # BB+RKHS
      M[,12]=f6b$yHat[w] # BB+3GK
      M[,13]=f7[w] # KNN
      M[,14]=predict(f8,gen[w,],150) # Boosting
      M[,15]=f18d # emRR
      M[,16]=f9b$hat[w] # Bag2
      M[,17]=f9c$hat[w] # Bag3
      M[,18]=f10$u[w] # REML
      M[,19]=f11[w] # GELA
      M[,20]=f12[w] # wGELA
      M[,21]=f13b$hat[w] # BagMEN
      M[,22]=predict(f14,gen[w,]) # EN
      M[,23]=pls16a # PLS
      M[,24]=f17a # varBB
      M[,25]=f17b # varBC
      M[,26]=f17c # SSVS
      M[,27]=f17d # EBL
      M[,28]=f17e # varMix
      M[,29]=f18a # emBA
      M[,30]=f18b # emBB
      M[,31]=f18c # emBC
      M[,32]=Y[w] # OBSERVARIONS
      return(M)
    }
    # Running Cross-validations
    b = list()
    for(i in 1:cvs){
      b[[i]] = folds(Seeds[i])
      cat('Done with',i,'of',cvs,'\n')
    }
    rm(list = list.files(pattern = 'ETA'))
    names(b) = paste('CV_',1:length(b),sep='')
    return(b)
  }
  
  CV.par=function(y,gen,k,Seeds,cores){
    # Check whether you have the packages
    check = function(mypkg) if(!is.element(mypkg,installed.packages()[,1])) install.packages(mypkg)
    check('BGLR')
    check('glmnet')
    check('gbm')
    check('rrBLUP')
    check('pls')
    check('kernlab')
    check('randomForest')
    check('snow')
    check('VIGoR')
    # Loading packages
    require(snow,quietly = T)
    # KNN function
    knn = function(y,E,k=20){d=function(d) mean(r[which(d<=(sort(d)[k]))]);
    w = which(is.na(y));D=E[,-w];r=y[-w];return(apply(D,1,d))}
    # Kernels
    E2 = as.matrix(dist(gen)^2)
    K=exp(-(E2/mean(E2)))
    K1=exp(-(E2/(quantile(E2,0.05))))
    K2=exp(-(E2/(0.2*quantile(E2,0.05))))
    K3=exp(-(E2/(5*quantile(E2,0.05))))
    eK = eigen(K,symmetric=T)
    eK1 = eigen(K1,symmetric=T)
    eK2 = eigen(K2,symmetric=T)
    eK3 = eigen(K3,symmetric=T)
    cat('DONE with EIGENDECOMPOSITIONS\n')
    Z = apply(gen,2,function(x) x-mean(x))
    G = tcrossprod(Z)
    G = G/mean(diag(G))
    Y = y
    N = nrow(gen)
    cvs = length(Seeds) 
    rm(Z,K1,K2,K3)
    # Cross-validation function
    folds = function(Seed,k,N,gen,G,eK,eK1,eK2,eK3,Y,y,E2,knn){
      # Loading packages
      require(BGLR,quietly = T)
      require(randomForest,quietly = T)
      require(glmnet,quietly = T)
      require(gbm,quietly = T)
      require(rrBLUP,quietly = T)
      require(pls,quietly = T)
      require(kernlab,quietly = T)
      require(VIGoR,quietly = T)
      # Begin folds
      set.seed(Seed)
      Nk = round(N/k)
      w=sample(1:N,Nk)
      y[w]=NA
      # RF
      f1=randomForest(x=gen[-w,],y=y[-w], mtry = min( floor(ncol(gen)/3), 1000 ))
      cat('RF\n')
      # Bayesian
      f2b = BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BayesB')),verbose=F)
      f2c = BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BayesC')),verbose=F)
      # RKHS
      f3 = BGLR(rmExistingFiles = F,y,ETA=list(list(V=eK$vectors,d=eK$values,model='RKHS')),verbose=F)
      f3b = BGLR(rmExistingFiles = F,y,ETA=list(list(V=eK1$vectors,d=eK1$values,model='RKHS'),
                                                list(V=eK2$vectors,d=eK2$values,model='RKHS'),
                                                list(V=eK3$vectors,d=eK3$values,model='RKHS')),verbose=F)
      cat('Alphabet and HS\n')
      # Lasso
      cv1 = cv.glmnet(x=gen[-w,],y=y[-w],alpha=1)
      lmb1 = cv1$lambda.min
      f4 = glmnet(x=gen[-w,],y=y[-w],lambda = lmb1,alpha = 1)
      f4b = BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BL')),verbose=F)
      # Ridge
      cv2 = cv.glmnet(x=gen[-w,],y=y[-w],alpha=1)
      lmb2 = cv2$lambda.min
      f5 = glmnet(x=gen[-w,],y=y[-w],lambda = lmb2,alpha = 0)
      f5b = BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BRR')),verbose=F)
      # BB+K
      f6 = BGLR(rmExistingFiles = F,y,ETA=list(list(V=eK$vectors,d=eK$values,model='RKHS'),
                                               list(X=gen,model='BayesB')),verbose=F)
      f6b = BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BayesB'),
                                                list(V=eK1$vectors,d=eK1$values,model='RKHS'),
                                                list(V=eK2$vectors,d=eK2$values,model='RKHS'),
                                                list(V=eK3$vectors,d=eK3$values,model='RKHS')),verbose=F)
      cat('done with L1 and L2\n')
      # kNN
      f7 = knn(y,E2,20)
      # Boosting
      f8=gbm::gbm.fit(x=gen[-w,],y=y[-w],distribution="gaussian",n.trees = 150)
      # Bag Bayes
      f9=wgr(y,gen,it=1800,bi=0,bag=.5,rp=T,pi=0.01,verb=T)#,eigK=eK,EigT = 1)
      f9b=wgr(y,gen,it=1800,bi=0,bag=.3,verb=T,rp=FALSE)
      f9c=wgr(y,gen,it=1800,bi=0,bag=.7,verb=T,rp=FALSE)
      # REML GBLUP
      f10 = mixed.solve(y,K=G)
      cat('done with Bags and GBLUP\n')
      # GELA
      f11 = GELA(y,x=gen,weighted = FALSE)
      f12 = GELA(y,gen,weighted = TRUE)
      cat('done with GELAs\n')
      # BagMAN
      f13b = ben(y=y,gen=gen,wpe = 25)
      # EN
      cv3 = cv.glmnet(x=gen[-w,],y=y[-w],alpha=0.5)
      lmb3 = cv3$lambda.min
      f14 = glmnet(x=gen[-w,],y=y[-w],lambda = lmb1,alpha = 0.5)
      # PLS'
      cat('now PLS\n')
      f16a=plsr(y[-w]~gen[-w,],validation='CV')
      ncomp=which.max(cor(f16a$validation$pred[,1,],y[-w]))
      pls16a = predict(f16a,gen[w,])[,1,ncomp]
      # RBF
      cat('now RBF\n')
      f2=predict(rvm(gen,y),gen)[,1]
      # varBayes
      f17a = gen[w,] %*% vigor(y,gen,'BayesB',hyperpara(gen,.5,'BayesB',.1),Function = 'tuning')$Beta
      f17b = gen[w,] %*% vigor(y,gen,'BayesC',hyperpara(gen,.5,'BayesC',.1),Function = 'tuning')$Beta
      f17c = gen[w,] %*% vigor(y,gen,'SSVS',hyperpara(gen,.5,'SSVS',.1),Function = 'tuning')$Beta
      f17d = gen[w,] %*% vigor(y,gen,'EBL',hyperpara(gen,.5,'EBL',.1),Function = 'tuning')$Beta
      f17e = gen[w,] %*% vigor(y,gen,'MIX',hyperpara(gen,.5,'MIX',.1),Function = 'tuning')$Beta
      # emBayes
      f18a = gen[w,] %*% emBA(y[-w],gen[-w,])$b
      f18b = gen[w,] %*% emBB(y[-w],gen[-w,])$b
      f18c = gen[w,] %*% emBC(y[-w],gen[-w,])$b
      f18d = gen[w,] %*% emRR(y[-w],gen[-w,])$b
      cat('Done with emBayes and varBayes\n')
      
      # Removing junk files and creating MODs
      file.remove(list.files(pattern='dat'))
      #mod = rep(NA,1+length(grep('^f.+',ls())))
      NamesMod = c('RandomForest','RBF','BayesB','BayesC','RKHS',
                   'RKHS(3GK)','LASSO','BL','RidgeReg',
                   'BRR','BayesB+RKHS','BayesB+3GK','kNN',
                   'Boost','emBRR','Bagging(30%)','Bagging(70%)',
                   'GBLUP-REML','GELA','wGELA','BEN','EN','PLS',
                   'varBayesB','varBayesC','SSVS','EBL','varMIX',
                   'emBayesA','emBayesB','emBayesC','OBSERVATION')
      M = matrix(NA,Nk,length(NamesMod))
      colnames(M) = NamesMod
      # VALUES
      M[,1]=predict(f1,gen[w,]) # RF
      M[,2]=f2[w] # RBF
      M[,3]=f2b$yHat[w] # BayesB
      M[,4]=f2c$yHat[w] # BayesC
      M[,5]=f3$yHat[w] # RKHS
      M[,6]=f3b$yHat[w] # RKHS - 3GK
      M[,7]=predict(f4,gen[w,]) # LASSO
      M[,8]=f4b$yHat[w] # BL
      M[,9]=predict(f5,gen[w,]) # Ridge
      M[,10]=f5b$yHat[w] # BRR
      M[,11]=f6$yHat[w] # BB+RKHS
      M[,12]=f6b$yHat[w] # BB+3GK
      M[,13]=f7[w] # KNN
      M[,14]=predict(f8,gen[w,],150) # Boosting
      M[,15]=f18d # emRR
      M[,16]=f9b$hat[w] # Bag2
      M[,17]=f9c$hat[w] # Bag3
      M[,18]=f10$u[w] # REML
      M[,19]=f11[w] # GELA
      M[,20]=f12[w] # wGELA
      M[,21]=f13b$hat[w] # BagMEN
      M[,22]=predict(f14,gen[w,]) # EN
      M[,23]=pls16a # PLS
      M[,24]=f17a # varBB
      M[,25]=f17b # varBC
      M[,26]=f17c # SSVS
      M[,27]=f17d # EBL
      M[,28]=f17e # varMix
      M[,29]=f18a # emBA
      M[,30]=f18b # emBB
      M[,31]=f18c # emBC
      M[,32]=Y[w] # OBSERVARIONS
      return(M)
    }
    # Create sock cluster
    cat('Creating Cluster with',cores,'cores\n')
    cl = makeSOCKcluster(cores)
    cat('Starting parallelization\n')
    b = clusterApply(cl, Seeds, fun=folds,k=k,N=N,gen=gen,G=G,eK=eK,eK1=eK1,eK2=eK2,eK3=eK3,Y=Y,y=y,E2=E2,knn=knn)
    cat('Closing cluster\n')
    stopCluster(cl)
    cat('All done\n')
    names(b) = paste('CV_',1:length(b),sep='')
    return(b)
  }
  
  if(is.null(cores)){
    H = CV.test(y=y,gen=gen,k=k,Seeds=Seeds)
  }else{
    H = CV.par(y=y,gen=gen,k=k,Seeds=Seeds,cores=cores)
  }
  
  return(H)
  
}

CV.check=function(cv){
  n = length(cv)
  m = ncol(cv$CV_1)
  dta = matrix(0,0,m)
  for(i in 1:n) dta = rbind(dta,cv[[i]])
  # functions
  MSE = function(A,B) mean((A-B)^2)
  TOPS = function(A,B,TOP=0.2) mean(which(B>quantile(B,1-TOP,na.rm = TRUE))%in%which(A>quantile(A,1-TOP,na.rm = TRUE)))
  BIAS = function(A,B) cov(A,B)/var(A)
  # summary
  PA = sort(cor(dta)[-m,m],decreasing = TRUE)
  Rank = sort(cor(dta,method = 'sp')[-m,m],decreasing = TRUE)
  MSPE = sort(apply(dta[,-m],2,MSE,B=dta[,m]))
  Top20 = sort(apply(dta[,-m],2,TOPS,B=dta[,m]),decreasing = TRUE)
  Bias = apply(dta[,-m],2,BIAS,B=dta[,m])
  Bias = Bias[order(as.matrix(dist(c(1,Bias)))[-1,1])]
  # Model choice
  o1 = o2 = o3 = o4 = o5 = 1:(m-1)
  names(o1) = names(PA); o1 = o1[order(names(o1))]
  names(o2) = names(PA); o2 = o2[order(names(o2))]
  names(o3) = names(PA); o3 = o3[order(names(o3))]
  names(o4) = names(PA); o4 = o4[order(names(o4))]
  names(o5) = names(PA); o5 = o5[order(names(o5))]
  O = o1+o2+o3+o4+o5
  ModelChoice = names(sort(O))[1:5]
  # output
  final = list('PA'=PA,'Rank'=Rank,
               'MSPE'=MSPE,'Bias'=Bias,
               'Top20'=Top20,'ModelChoice'= ModelChoice)
  return(final)
  
}
