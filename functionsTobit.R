
# ----------------------------- EXACT iid SAMPLER ------------------------------
rSUNpost = function(X1,y1,X0,y0,om2,zT,nSample,seed=NULL){
  
  print(" - iid Sampler from the exact SUN posterior - ")
  
  nSample = round(nSample/10)*10
  
  ### algebraic precomputations
  print(" Precomputations ")
  
  if(!is.null(seed)){set.seed(seed)}
  
  n0 = dim(X0)[1]
  n1 = dim(X1)[1]
  p = dim(X0)[2]
  
  if(p<n1){
    Omega_post = solve(diag(1./om2,p,p) + crossprod(X1))
  }else{
    Omega_post = diag(om2,p,p) - om2*crossprod(X1,solve(diag(1.,n1,n1)+om2*tcrossprod(X1))%*%X1)*om2
  }
  Omega_post = 0.5*(Omega_post+t(Omega_post))
  xi_post = Omega_post%*%(crossprod(X1,as.matrix(y1)) + c(-zT/om2,rep(0,p-1)))
  
  signX0 = X0
  signX0[y0==0,] = -X0[y0==0,]
  
  Gamma_post_unnormalized = diag(1.,n0,n0)+signX0%*%tcrossprod(Omega_post,signX0)
  inv_Gamma_post_unnormalized = solve(Gamma_post_unnormalized)
  s = sqrt(Gamma_post_unnormalized[cbind(1:n0,1:n0)])
  s_1 = 1/s
  gamma_post = s_1*(signX0%*%xi_post)
  Gamma_post = s_1*t(s_1*Gamma_post_unnormalized)
  
  ### compute multiplicative coefficients for the truncated multivariate normal component
  X0Omega = signX0%*%Omega_post
  invGuX0Omega = inv_Gamma_post_unnormalized%*%X0Omega
  
  ### compute multiplicative coefficients for the multivariate normal component
  V = Omega_post-crossprod(X0Omega,invGuX0Omega)
  V = 0.5*(V+t(V))
  L = t(chol(V))
  
  ### sample the multivariate normal component
  print(" Sampling Multivariate Normal ")
  sampleMultNorm = matrix(rnorm(nSample*p),p,nSample)
  
  sampleMN = L%*%sampleMultNorm
  
  ### sample the multivariate truncated normal component
  print(" Sampling Truncated Normal ")
  sampleTruncNorm = matrix(,nrow=n0,ncol=0)
  for(i in c(1:10)){
    print(sprintf(fmt = "%10s%3s%2s", "[",(i-1)*10,"%]"))
    sampleTruncNorm = cbind(sampleTruncNorm,t(rtmvnorm(n = nSample/10, mu = rep(0,n0), sigma = Gamma_post, lb = -gamma_post, u = rep(Inf,n0))) )
  }
  print(sprintf(fmt = "%10s%3s%2s", "[",i*10,"%]"))
  
  sampleTN = crossprod(s*invGuX0Omega,sampleTruncNorm)
  
  ### combine the multivariate normal and truncated normal components
  print(" Combining Results ")
  sampleSUN = matrix(xi_post,p,nSample) + sampleMN + sampleTN
  
  results = list(sampleSUN = sampleSUN)
  
  return(results)
}



# ----------------------------------- MF-VB ------------------------------------
getParamsMF = function(X1,y1,X0,y0,om2,zT,tolerance=1e-3,
                       maxIter=1e3,fullVar=FALSE,nPrint=1000){
  
  n1 = dim(X1)[1]
  n0 = dim(X0)[1]
  p = dim(X0)[2]
  
  ### Pre-Computations
  
  r0 = crossprod(X1,y1) + c(-zT/om2,rep(0,p-1))
  if(p<n1){
    Omega0 = solve(diag(1/om2,p,p) + crossprod(X1))
    beta0 = Omega0%*%r0
  }else{
    invXXt = solve(diag(1.,n1,n1)+om2*tcrossprod(X1))
    Omega0 = diag(om2,p,p) - om2*t(X1)%*%invXXt%*%X1*om2
    beta0 = om2*(r0 - crossprod(X1,invXXt%*%(X1%*%r0))*om2)
  }
  # Omega0 = 0.5*(Omega0+t(Omega0))
  
  if(p<n0) {
    invOmega0 = solve(Omega0)
    V = solve(crossprod(X0)+invOmega0)
    diagV = diag(V)
    VXt = tcrossprod(V,X0)
    XVinvOmega0 = crossprod(VXt,invOmega0)
  } else{
    Omega0Xt = tcrossprod(Omega0,X0)
    Lambda = solve(diag(1, nrow = n0, ncol = n0)+X0%*%Omega0Xt)
    VXt = Omega0Xt%*%Lambda
    XVinvOmega0 = Lambda%*%X0
    diagV = diag(Omega0) - rowSums(Omega0Xt * VXt)
    H = XVinvOmega0%*%VXt # needed for ELBO
  }

  VinvOmega0beta0 = beta0 - VXt%*%(X0%*%beta0)
  XVinvOmega0beta0 = XVinvOmega0%*%beta0
  rELBO = XVinvOmega0%*%VinvOmega0beta0 + XVinvOmega0beta0
  
  if(p<n0){
    signVXt = VXt
    signVXt[,y0==0] = -signVXt[,y0==0]
  } else {
    XVXt = X0%*%VXt
    signXVXt = XVXt
    signXVXt[y0==0,] = -XVXt[y0==0,]
  } 
  
  signXVinvOmega0beta0 = XVinvOmega0beta0
  signXVinvOmega0beta0[y0==0,] = -signXVinvOmega0beta0[y0==0,]
  
  ### Initialization
  
  mu = matrix(0,n0,1)
  meanZ = mu + (2*y0-1)*exp(dnorm(mu, log = T) - pnorm((2*y0-1)*mu, log.p = T))
  diff = 1
  elbo = -Inf
  nIter = 0
  
  ### Iterations
  
  if(p<n0){
    
    while(diff > tolerance & nIter < maxIter) {
    
      mu = crossprod(VXt,crossprod(X0,meanZ)) + XVinvOmega0beta0
      meanZ = mu + (2*y0-1)*exp(dnorm(mu, log = T) - pnorm((2*y0-1)*mu, log.p = T))
      
      elboOld = elbo
      
      elbo = -sum(crossprod(XVinvOmega0,meanZ) * (VXt%*%meanZ))/2 + sum(meanZ*rELBO) +
        sum(pnorm(crossprod(signVXt,crossprod(X0,meanZ))+signXVinvOmega0beta0, log.p = T))
      
      diff = abs(elbo-elboOld)
      
      nIter = nIter+1
      if(nIter %% nPrint == 0) {print(paste0("iteration: ", nIter, ", ELBO: ", elbo))}
    }
    
  } else {
    
    while(diff > tolerance & nIter < maxIter) {
      
      mu = XVXt%*%meanZ + XVinvOmega0beta0
      meanZ = mu + (2*y0-1)*exp(dnorm(mu, log = T) - pnorm((2*y0-1)*mu, log.p = T))
      
      elboOld = elbo
      
      elbo = -sum(meanZ *(H%*%meanZ))/2 + t(meanZ)%*%rELBO +
        sum(pnorm(signXVXt%*%meanZ+signXVinvOmega0beta0, log.p = T))
      
      diff = abs(elbo-elboOld)
      
      nIter = nIter+1
      if(nIter %% nPrint == 0) {print(paste0("iteration: ", nIter, ", ELBO: ", elbo))}
    }
    
  }
  
  ### Posterior Approximate Moments
  
  meanBeta = VXt%*%meanZ + VinvOmega0beta0
  
  results = list(meanBeta = meanBeta, diagV = diagV, nIter = nIter)
  
  if(fullVar==TRUE){
    if(p>=n0) {
      V = Omega0 - tcrossprod(Omega0Xt,VXt)
    }
    results = append(list(V=V),results)
  }
  
  return(results)
}




# ---------------------------------- PFM-VB ------------------------------------
getParamsPFM = function(X1,y1,X0,y0,om2,zT,tolerance=1e-3,maxIter=1e3,
                        fullVar=FALSE,predictive=FALSE,nPrint=1000){
  
  n1 = dim(X1)[1]
  n0 = dim(X0)[1]
  p = dim(X0)[2]
  
  ### Pre-Computations
  
  r0 = crossprod(X1,y1) + c(-zT/om2,rep(0,p-1))
  if(p<n1){
    Omega0 = solve(diag(1/om2,p,p) + crossprod(X1)) 
    beta0 = Omega0%*%r0
  }else{
    invXXt = solve(diag(1.,n1,n1)+om2*tcrossprod(X1))
    Omega0 = diag(om2,p,p) - om2*t(X1)%*%invXXt%*%X1*om2
    beta0 = om2*(r0 - crossprod(X1,invXXt%*%(X1%*%r0))*om2)
  }
  # Omega0 = 0.5*(Omega0+t(Omega0))
    
  if(p<n0) {
    invOmega0 = solve(Omega0)
    V = solve(crossprod(X0)+invOmega0)
    H = X0%*%tcrossprod(V,X0)
    Lambda = -H    # more efficient than Lambda = diag(1,nrow=n0,ncol=n0) - H
    Lambda[cbind(1:n0,1:n0)] = 1+Lambda[cbind(1:n0,1:n0)]
  } else{
    Omega0Xt = tcrossprod(Omega0,X0)
    Lambda = solve(diag(1,nrow=n0,ncol=n0)+X0%*%Omega0Xt)
    H = diag(1,nrow=n0,ncol=n0)-Lambda
  }
  
  Xbeta0 = X0%*%beta0 
  LambdaXbeta0 = Lambda%*%Xbeta0
  
  sigma2 = as.vector(1/(1-diag(H)))
  sigma = as.vector(sqrt(sigma2))
  
  A = sigma2*H
  A[cbind(1:n0,1:n0)] = 0

  ### Initialization
  
  mu = double(n0)
  LogPhi = double(n0)
  
  musiRatio = (2*y0-1)*mu/sigma
  phiPhiRatio = exp(dnorm(musiRatio,log = T)-pnorm(musiRatio,log.p = T))
  meanZ = mu + (2*y0-1)*sigma*phiPhiRatio
  
  if(p<n0){
    XV = X0%*%V
    alpha = crossprod(meanZ,X0) - meanZ[n0]*X0[n0,]
    iMinus = c(n0,1:(n0-1))
  }
  
  elbo = -Inf
  diff = 1
  nIter = 0
  
  ### Iterations
  
  if(p<n0){
    
    while(diff > tolerance & nIter < maxIter) {
      
      elboOld = elbo
    
      for(i in 1:n0) {
        
        alpha = alpha + X0[iMinus[i],]*meanZ[iMinus[i]] - X0[i,]*meanZ[i]
        mu[i] = sigma2[i]*LambdaXbeta0[i] + sigma2[i]*alpha%*%XV[i,]
        
        musiRatio = (2*y[i]-1)*mu[i]/sigma[i]
        LogPhi[i] = pnorm(musiRatio, log.p = T)
        phiPhiRatio = exp(dnorm(musiRatio, log = T) - LogPhi[i])
        meanZ[i] = mu[i] + (2*y[i]-1)*sigma[i]*phiPhiRatio
        
      }
      
      XtmeanZ = crossprod(X0,meanZ)
      elbo = -sum(meanZ^2)/2 + sum(XtmeanZ*(V%*%XtmeanZ))/2 +
        sum(((meanZ-mu)^2)/sigma2)/2 + sum(LogPhi) + sum(meanZ*LambdaXbeta0)
    
      diff = abs(elbo-elboOld)
      nIter = nIter+1
    
      if(nIter%%nPrint==0) {print(paste0("iteration: ", nIter, ", ELBO: ", elbo))}
    }
  
  } else {
    
    while(diff > tolerance & nIter < maxIter) {
      
      elboOld = elbo
      
      for(i in 1:n0) {
        
        mu[i] = sigma2[i]*LambdaXbeta0[i] + A[i,]%*%meanZ
        
        musiRatio = (2*y[i]-1)*mu[i]/sigma[i]
        LogPhi[i] = pnorm(musiRatio, log.p = T)
        phiPhiRatio = exp(dnorm(musiRatio, log = T) - LogPhi[i])
        meanZ[i] = mu[i] + (2*y[i]-1)*sigma[i]*phiPhiRatio
      }
      
      elbo = -crossprod(meanZ,Lambda)%*%meanZ/2 + sum(((meanZ-mu)^2)/sigma2)/2 +
        sum(LogPhi) + sum(meanZ*LambdaXbeta0)
    
      diff = abs(elbo-elboOld)
      nIter = nIter+1
      
      if(nIter%%nPrint==0) {print(paste0("iteration: ", nIter, ", ELBO: ", elbo))}
    }
    
  }
  
  mu = Xbeta0+A%*%(meanZ-Xbeta0)
  
  results = list(mu = mu, sigma2 = sigma2, nIter = nIter)
  
  ### Posterior Approximate Moments
  
  if(p<n0) {
    diagV = diag(V) # V already computed
    VXt = tcrossprod(V,X0)
    VinvOmega0 = V%*%invOmega0
  } else{
    VXt = Omega0Xt%*%Lambda
    VinvOmega0 = diag(1,nrow=p,ncol=p) - VXt%*%X0
    diagV = diag(Omega0) - rowSums(Omega0Xt * VXt)

  }
  
  musiRatio = (2*y0-1)*mu/sigma
  phiPhiRatio = exp(dnorm(musiRatio, log = T) - pnorm(musiRatio, log.p = T))
  
  meanZ = mu + (2*y0-1)*sigma*phiPhiRatio
  postVarZ = as.double(sigma2*(1-musiRatio*phiPhiRatio - phiPhiRatio^2))
  
  W = rowSums(VXt*t(postVarZ*t(VXt)))
  
  VinvOmega0beta0 = VinvOmega0%*%beta0
  
  meanBeta = VXt%*%meanZ + VinvOmega0beta0
  varBeta = diagV + W
  
  moments_PFM = list(meanBeta=meanBeta,varBeta=matrix(varBeta,ncol = 1))
  
  if(fullVar==TRUE){
    if(p>=n0) {
      V = Omega0 + tcrossprod(VXt,Omega0Xt)
    }
    varBeta = V + VXt%*%(postVarZ*t(VXt))
    moments_PFM$varBeta=varBeta
  }
  
  results = c(results,postMoments=moments_PFM) 
  
  if(predictive==TRUE){
    for_preditive_locations = list(VXt=VXt,VinvOmega0beta0=VinvOmega0beta0)
    results = c(results,postPredictive=for_preditive_locations)
  }
  
  return(results)
}




# ------------------------------- EFFICIENT EP ---------------------------------
zeta1 = function(x){exp(dnorm(x, log = T) - pnorm(x,log.p = T))}
zeta2 = function(x,z1){-z1^2-x*z1}

getParamsEP = function(X1,y1,X0,y0,om2,zT,tolerance=1e-3,maxIter=1e3,
                       fullVar=FALSE,predictive=FALSE,nPrint=1000){
  
  n1 = dim(X1)[1]
  n0 = dim(X0)[1]
  p = dim(X0)[2]
  
  ### Pre-Computations
  
  r = crossprod(X1,y1) + c(-zT/om2,rep(0,p-1))
  
  if(p<n1){
    Omega0 = solve(diag(1./om2,p,p) + crossprod(X1))
    beta0 = Omega0%*%r
  }else{
    Lambda = solve(diag(1.,n1,n1)+om2*tcrossprod(X1))
    Omega0 = diag(om2,p,p) - om2*crossprod(X1,Lambda%*%X1)*om2
    beta0 = om2*(r - crossprod(X1,Lambda%*%(X1%*%r))*om2)
  }
  # Omega0 = 0.5*(Omega0+t(Omega0))
  
  if(p<n0){
    Omega = Omega0
    if(p<n1){
      Q = solve(Omega)
    }else{
      Q = diag(1./om2,p,p) + crossprod(X1)
    }
    logDetQ0 = determinant(Q, logarithm = TRUE)
  }else{
    if(p<n1){
      U = tcrossprod(Omega0,X0)
    } else {
      U = om2*t(X0) - om2*crossprod(X1,Lambda%*%tcrossprod(X1,X0))*om2
    }
    U0 = U
  }
  
  logZ0 = 0.5*sum(r*beta0)
  logZ = double(length = n0)
  
  ### Initialization
  
  k = double(length = n0)
  m = double(length = n0)
  
  diff = 1
  nIter = 0
  
  ### Iterations
  
  if(p<n0){
      
    while(diff > tolerance && nIter < maxIter){
        
      diff = 0.
      count = 0
      logZ = 0.
      
      for(i in c(1:n0)){
        
        xi = as.matrix(X0[i,])
        
        r_i = r - m[i]*xi
        Q_i = Q - k[i]*tcrossprod(xi)
        
        Oxi = Omega%*%xi
        
        Oi = Omega + tcrossprod(Oxi)*k[i]/as.double(1.-k[i]*crossprod(xi,Oxi))
        
        Oixi = Oi%*%xi
        xiOixi = as.double(crossprod(xi,Oixi))
        
        if(xiOixi>0){
          
          r_iOixi = as.double(crossprod(r_i,Oixi))
          
          s = (2.*y0[i]-1.)/sqrt(1.+xiOixi)
          tau = s*r_iOixi
          
          z1 = zeta1(tau)
          z2 = zeta2(tau,z1)
          
          kNew = - z2/(1.+xiOixi+z2*xiOixi)
          mNew = s*z1 + kNew*r_iOixi + kNew*s*z1*xiOixi
          
          maxDiff = max(abs(c(kNew - k[i], mNew - m[i])))
          if(maxDiff>diff){diff = maxDiff}
          
          k[i] = kNew
          m[i] = mNew
          
          logZ[i] = pnorm(tau,log.p = T) + 0.5*log(1.+k[i]*xiOixi) +
            0.5*((r_iOixi)^2)/xiOixi - 0.5*((m[i] + r_iOixi/xiOixi)^2)*xiOixi/(1.+k[i]*xiOixi)
          
          r = r_i + m[i]*xi
          Q = Q_i + k[i]*tcrossprod(xi)
          
          Omega = Oi - tcrossprod(Oixi)*k[i]/(1.+k[i]*xiOixi)
        }else{
          count = count+1
          print(paste0(count," units skipped"))
        }
        
      }
    
      nIter = nIter + 1
      if(nIter %% nPrint == 0) {print(paste0("iteration ",nIter))}
    }
    
  } else {
    
    while(diff > tolerance && nIter < maxIter){
      
      diff = 0.
      count = 0
      logZ = 0.
      
      for(i in c(1:n0)){
        
        u = U[,i]
        xi = X0[i,]
        xTu = sum(xi*u)
        
        d = 1-k[i]*xTu
        w = u/d
        xTw = xTu/d
        
        if(xTw>0){
          
          r_iTw = sum(r*w) - m[i]*xTw
          
          s = (2*y0[i]-1)/sqrt(1+xTw)
          tau = s*r_iTw
          
          z1 = zeta1(tau)
          z2 = zeta2(tau,z1)
          
          kNew = as.double(-z2/(1 + xTw + z2*xTw))
          mNew = as.double(z1*s + kNew*r_iTw + kNew*z1*s*xTw)
          
          r = r + (mNew - m[i])*xi
          
          c = (k[i]-kNew)/(1+(kNew-k[i])*xTu)
          U = U + (u*c)%*%crossprod(xi,U)
          
          maxDiff = max(abs(c(kNew - k[i], mNew - m[i])))
          if(maxDiff>diff){diff = maxDiff}
          
          k[i] = kNew
          m[i] = mNew
          
          logZ[i] = pnorm(tau,log.p = T) + 0.5*log(1.+k[i]*xTw) +
            0.5*((r_iTw)^2)/xTw - 0.5*((m[i]+r_iTw/xTw)^2)*xTw/(1.+k[i]*xTw)
          
        }else{
          count = count+1
          print(paste0(count," units skipped"))
        }
      }
      
      nIter = nIter + 1
      if(nIter %% nPrint == 0) {print(paste0("iteration ",nIter))}
    }
  
  }
  
  ### Posterior Approximate Moments
  
  if(p<n0){
    
    meanBeta = Omega%*%r
    diagOmega = diag(Omega)
    
    logDetQ = determinant(Q, logarithm = TRUE)
    logML = sum(logZ) + 0.5*sum(r*meanBeta) - logZ0 + 
      0.5*logDetQ0$modulus[1] - 0.5*logDetQ$modulus[1] 
    
  }else{
    
    diagOmega = diag(Omega0) - rowSums(U*t(k*t(U0)))
    meanBeta = beta0 + U0%*%m - U%*%(k*crossprod(U0,r))
    
    logDet = determinant(diag(1,n0,n0) + k*X0%*%U0, logarithm = TRUE)
    logML = sum(logZ) + 0.5*sum(r*meanBeta) - logZ0 - 0.5*logDet$modulus[1]
  }
  
  results = list(meanBeta = meanBeta, diagOmega = diagOmega, logML = logML, 
                 nIter = nIter, kEP = k, mEP = m)
  
  if(fullVar==TRUE){
    if(p>=n0) {
      Omega = Omega0 - U%*%(k*t(U0))
    }
    results = append(list(Omega=Omega),results)
  }
  
  if(predictive==TRUE){
    if(p>=n0) {
      results = c(results,postPredictive=list(U=U,U0=U0))
    }
  }
  
  return(results)
}



# -------------------------- GENERATE SIMULATED DATA ---------------------------
generateData = function(n,p,k,nTest=200,beta=NULL,seed=123,input_dir="dataset"){
  
  set.seed(seed)
  
  # True Parameters
  if(is.null(beta)){
    beta = matrix(runif(p,-5,5), ncol = 1)
  }
  
  # Design Matrix
  X = matrix(rnorm((p-1)*(n+nTest),0,1), nrow = n+nTest, ncol = p-1)
  X <- apply(X,2,function(y) y - mean(y))
  X <- apply(X,2,arm::rescale)
  X = cbind(rep(1,n+nTest),X)
  
  # Latent Utilities
  z = X%*%beta + matrix(rnorm(n+nTest),n+nTest,1)
  
  # Censoring Threshold
  zT = as.double(quantile(z, c(k/100)))
  
  # Observed Responses
  y = pmax(rep(0,n+nTest),z-zT)
  
  # Censored observations
  y0 = y[y==0]
  X0 = X[y==0,]
  
  # Uncensored observations
  y1 = y[y>0]
  X1 = X[y>0,]
  
  # Divide the data in (proportioned) train and test observations
  testCens = ceiling(k*nTest/100)
  
  # Out of sample
  Xtest = rbind(X0[1:testCens,],X1[1:(nTest-testCens),])
  yTest = c(y0[1:testCens],y1[1:(nTest-testCens)])
  
  # In sample
  X = rbind(X0[-(1:testCens),],X1[-(1:(nTest-testCens)),])
  y = c(y0[-(1:testCens)],y1[-(1:(nTest-testCens))])
  
  save(X,y,nTest,Xtest,yTest,zT,
       file=paste(input_dir,"/simulatedData_n",n,"_p",p,"_k",k,".RData",sep=""))
}

