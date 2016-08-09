
wgr = function(y,gen,it=1500,bi=500,th=1,bag=1,rp=FALSE,iv=FALSE,pi=0,df=5,R2=0.5,eigK=NULL,rankK=0.25,verb=FALSE){
  
  anyNA = function(x) any(is.na(x))
  if(anyNA(gen)){stop('No missing values allowed in Z')}
  if(anyNA(pi>=1)){stop('Pi must be lower than 1')}
  gen0 = gen
  
  # Polygenic term
  if(!is.null(eigK)){
   pk = round(rankK*length(eigK$values))
   U0 = U = eigK$vectors[,1:pk]
   V = eigK$values[1:pk]
   H = h = rep(0,pk)
   Vk = rep(1,pk)
   xxK = rep(bag,pk)
  }  
  
  # Remove missing y's
  if(anyNA(y)){
    mis = which(is.na(y))
    y = y[-mis]
    gen = gen[-mis,]
    if(!is.null(eigK)) U = U[-mis,]
  }
  
  # MCMC settings
  post = seq(bi,it,th)
  mc = length(post)
  
  # Preparing inputs
  X = gen
  n = nrow(gen)
  p = ncol(gen)
  xx = apply(X,2,crossprod)*bag
  b = g = rep(0,p)
  d = rep(1,p)
  mu = mean(y - X%*%b)
  e = y - X%*%g - mu
  Va = 0.0001
  Vb = rep(Va,p)
  Ve = 1
  md = pi
  
  # Priors
  df_prior = df 
  MSx = sum(apply(X,2,var,na.rm=T))
  S_prior = R2*var(y,na.rm=T)*(df_prior+2)/MSx/(1-pi)
  shape_prior  = 1.1
  rate_prior = (shape_prior-1)/S_prior
  S_conj = rgamma(1,p*df_prior/2+shape_prior,sum(1/Vb)/2+rate_prior)
  if(!is.null(eigK)) Sk_prior = R2*var(y,na.rm=T)*(df_prior+2)/pk
  
  # Storing Posterior
  B0 = VA = VE = VP = S = 0
  VB = G = D = B = rep(0,p)
  
  #RUN
  if(verb) pb = txtProgressBar(style = 3)
  
  for(i in 1:it){
    
    # Resampling
    if(bag!=1) Use = sort(sample(n,n*bag,rp))
    
    # Update polygenic term and regression coefficients
    if(!is.null(eigK)){
      
      Lk = Ve/(V*Vk)
      if(bag!=1){ update = KMUP(U[Use,],h,xxK,e[Use],Lk,pk,Ve,0)
      }else{ update = KMUP(U,h,xxK,e,Lk,pk,Ve,0)}
      h = update[[1]]
      e = update[[3]]
      
      if(pi>0){PI = rbeta(1,10*pi+md+1,10*(1-pi)-md+1)
      }else{PI=0}  
      L = Ve/Vb
      if(bag!=1){update = KMUP(X[Use,],g,xx,e,L,p,Ve,PI)
      }else{ update = KMUP(X,g,xx,e,L,p,Ve,PI)}
      #if(bag!=1){update = KMUP(X[Use,],b,xx,e,L,p,Ve,PI)
      #}else{ update = KMUP(X,b,xx,e,L,p,Ve,PI)}
      b = update[[1]]
      d = update[[2]]; if(pi>0) d[is.nan(d)] = 1
      e = update[[3]]
      g = b*d
      md = mean(d)
      
    }else{
    
      # Update regression coefficients without polygene
      if(pi>0){PI = rbeta(1,25*pi+md+1,25*(1-pi)-md+1)}else{PI=0}  
      L = Ve/Vb
      #if(bag!=1){update = KMUP(X[Use,],b,xx,e[Use],L,p,Ve,PI)
      #}else{ update = KMUP(X,b,xx,e,L,p,Ve,PI)}
      if(bag!=1){update = KMUP(X[Use,],g,xx,e[Use],L,p,Ve,PI)
      }else{ update = KMUP(X,g,xx,e,L,p,Ve,PI)}
      b = update[[1]]
      d = update[[2]]; if(pi>0) d[is.nan(d)] = 1
      e = update[[3]]
      g = b*d
      md = mean(d)
      
    }
            
    # Update genetic variance
    if(iv){    
      Vb = (S_conj+b^2)/rchisq(p,df_prior + 1)
      S_conj = rgamma(1,p*df_prior/2+shape_prior,sum(1/Vb)/2+rate_prior)
    }else{
      Va = (sum(b^2)+S_prior)/rchisq(1,df_prior + p)
      Vb = rep(Va,p)
    }
    if(!is.null(eigK)){
      Vp = (sum(h^2/V)+Sk_prior)/rchisq(1,df_prior + pk)
      Vk = rep(Vp,pk)
    }
    
    
    # Residual variance
    Ve = (sum(e*e)+S_prior)/rchisq(1,n*bag + df_prior)
        
    # Intercept
    if(!is.null(eigK)){e = y-X%*%g-U%*%h}else{e = y-X%*%g}
      mu = rnorm(1,mean(e),(Ve+1E-5)/n)
      if(is.na(mu)||is.nan(mu)) mu = mean(y,na.rm=T)
    e = e-mu

    # Save posterior
    if(i%in%post){
      B0 = B0+mu
      B = B+b
      D = D+d
      G = G+g
      VE = VE+Ve
      if(iv){VB = VB+Vb}else{VA = VA+Va}
      if(!is.null(eigK)){H = H+h; VP = VP+Vp}
    }
    
    if(verb) setTxtProgressBar(pb, i/it)
  }
  
  if(verb) close(pb)
  
  # Posterior mean
  B0 = B0/mc
  B = B/mc
  D = D/mc
  G = G/mc
  VE = VE/mc
  if(iv){VB = VB/mc}else{VA = VA/mc}
  if(!is.null(eigK)){H = H/mc; VP = VP/mc}
  
  # Fitted values
  if(!is.null(eigK)){
    poly = U0 %*% H
    HAT = B0 + gen0 %*% G + poly
  }else{
    HAT = B0 + gen0 %*% G
  }
  
  # Output
  if(!is.null(eigK)){
    if(iv){
      final = list('mu'=B0,'b'=B,'Vb'=VB,'g'=G,'d'=D,'Ve'=VE,'hat'=HAT,'u'=poly,'Vk'=VP)
    }else{
      final = list('mu'=B0,'b'=B,'Vb'=VA,'g'=G,'d'=D,'Ve'=VE,'hat'=HAT,'u'=poly,'Vk'=VP)
    }
  }else{
    if(iv){
      final = list('mu'=B0,'b'=B,'Vb'=VB,'g'=G,'d'=D,'Ve'=VE,'hat'=HAT)
    }else{
      final = list('mu'=B0,'b'=B,'Vb'=VA,'g'=G,'d'=D,'Ve'=VE,'hat'=HAT)
    }
  }
  

  
  return(final)
}

ben = function(y,gen,it=750,bi=250,bag=0.5,alpha=0.5,wpe=10){
  X=gen;
  rm(gen);
  
  # Function to update beta
  upB = function(e,mu,X,b,l,a,xx,p,E2,X2,bag,pi,wpe,O){
    xx = xx*bag
    a_new = rbeta(1,10*pi,10*(1-pi))
    b1 = b2 = rep(NA,p)
    e1 = e2 = e
    L1_1 = (1-a)*l/(2*xx)
    L2_1 = (a)*(xx/(xx*l))
    L1_2 = (1-a_new)*l/(2*xx)
    L2_2 = (1-a_new)*(xx/(xx*l))
    # Update regression coefficients
    for(j in O){
      # Old alpha
      xy1 = (crossprod(e1,X[,j])+b[j]*xx[j])/xx[j]
      s1 = sign(xy1)
      beta = abs(xy1)-L1_1[j]
      ben = s1*beta*L2_1[j]
      b1[j] = ifelse(beta>0,ben,0)
      e1 = e1 - X[,j]*(b1[j]-b[j])
      # New alpha
      xy2 = (crossprod(e2,X[,j])+b[j]*xx[j])/xx[j]
      s2 = sign(xy2)
      beta = abs(xy2)-L1_2[j]
      ben = s2*beta*L2_2[j]
      b2[j] = ifelse(beta>0,ben,0)
      e2 = e2 - X[,j]*(b2[j]-b[j])
    }
    # Loss function
    SSPE_1 = sum(as.vector(tcrossprod(b1,X2)-E2)^2)+0.00001
    SSPE_2 = sum(as.vector(tcrossprod(b2,X2)-E2)^2)+0.00001
    LOSS1 = wpe*SSPE_1+crossprod(e1)+l*(0.5*crossprod(b1)*(1-a)+sum(abs(b1))*a)
    LOSS2 = wpe*SSPE_2+crossprod(e2)+l*(0.5*crossprod(b2)*(1-a_new)+sum(abs(b2))*a_new)
    LR = LOSS2/(LOSS1+LOSS2)
    # METROPOLIS
    if(LR<runif(1)){
      P=list('b'=b2,'a'=a_new,'e'=e2,'oob'=SSPE_2)
    }else{
      P=list('b'=b1,'a'=a,'e'=e1,'oob'=SSPE_2)}
    return(P)
  }
  # Missing
  if(anyNA(y)){
    mis = which(is.na(y))
    Y = y
    XX = X[mis,]
    y = y[-mis]
    X = X[-mis,]
    MISS = TRUE
  }else{
    MISS = FALSE
  }
  # Data
  xx = apply(X,2,function(x)crossprod(x))
  b0 = crossprod(X,y)[,1]/xx
  #SP = 0.5*(4+2)*var(y,na.rm=T)/sum(apply(X,2,var,na.rm=T))
  O = order(b0^2,decreasing = TRUE)
  n = nrow(X)
  p = ncol(X)
  bn = round(n*bag)
  MC = it-bi
  # Parameters
  mu = mean(y)
  e = y-mu
  b = rep(0,p)
  a = 1
  l = 1
  # Store posterior
  B = rep(0,p)
  MU = A = L = SSPE = 0
  # Loop
  pb = txtProgressBar(style = 3)
  for(i in 1:it){
    
      # Bagging
      s = sort(sample(1:n,n-bn,replace=FALSE))
      # UPDATE
      UP = upB(e[-s],mu,X[-s,],b,l,a,xx*bag,p,e[s],
               X[s,],bag,alpha,wpe,O=O)
      b = UP[[1]]
      a = UP[[2]]
      e = UP[[3]]
    
    mu = mu + mean(e)
    S_prior = runif(1)#*SP
    df_prior = rpois(1,4)
    ve = crossprod(e+S_prior)/(bn+df_prior)
    vb = crossprod(b+S_prior)/(p+df_prior)
    l = ve/vb
    e = as.vector(y-(mu+tcrossprod(b,X)))
    # STORE
    if(it>bi){
      B = B+b
      MU = MU+mu
      A = A+a
      L = L+l
      if(bag<1) SSPE = SSPE+UP$oob/(n-bn)
    }
    setTxtProgressBar(pb, i/it)
  }
  close(pb)
  # Posterior
  Bhat = B/MC
  MUhat = MU/MC
  Ahat = A/MC
  Lhat = L/MC
  MSPEout = mean(SSPE)/MC
  # Prediction
  if(MISS){
    Yhat = Y
    Yhat[-mis] = as.vector(mu+tcrossprod(Bhat,X))
    Yhat[mis] = as.vector(mu+tcrossprod(Bhat,XX))
  }else{
    Yhat = as.vector(mu+tcrossprod(Bhat,X))
  }
  # OUTPUT
  LIST = list('hat'=Yhat,'coef'=Bhat,'b0'=MUhat,
              'alp'=Ahat,'lmb'=Lhat,'MSPEoob'=MSPEout)
  
}
