
if(sys.nframe()==0L){rm(list = ls())}
source("./covglasso_ada.R")

LcUc<-function(lower, upper, rr, c){
  nnc = length(c)# number of censored variables
  Lc <- Uc<- rep(0, nnc)
  j <- 0
  for(i in c){
    j = j+1
    if(rr[i]==1){
      Lc[j] = upper[i]
      Uc[j] = Inf
    }
    else if(rr[i]==-1){
      Lc[j] = -Inf
      Uc[j] = lower[i]
    }
  }
  list(lower=Lc, upper=Uc)
}
cond.mean.Precision <-function(mu, H, o, c, cX){
  
  if(length(o)>0){
    
    if(length(c) == 1){
      muc = mu[c] - solve(H[c,c])%*% H[c,o]%*%matrix(cX[o] - mu[o], ncol=1)
    }
    else{
      muc = mu[c] - solve(H[c,c], H[c,o])%*%matrix(cX[o] - mu[o], ncol=1)
    }
    
    Hcc = H[c,c]
    
    
    if(length(c)>1){
      Hcc[abs(Hcc)<1e-7]<-0
      Hcc = as(Hcc, Class="dgCMatrix")# sparseMatrix(i = ep, j = ep, x=, dims=c(,))
    }
  }
  else{
    if(length(c) == length(mu)){
      muc = mu
      Hcc = H
      Hcc[abs(Hcc)<1e-7]<-0
      Hcc = as(Hcc, Class="dgCMatrix")# sparseMatrix(i = ep, j = ep, x=, dims=c(,))
    }
    else
      stop("length(o) = 0, but length(c) != length(mu)")
    
  }
  
  return(list(muc=muc, Hcc=Hcc))
}

compute.loglik<-function(n, S, wi){
  
  -0.5*n*nrow(S)*log(2*pi) + 0.5*n*log(det(wi)) - 0.5*n*sum(diag(S%*%wi))
}

mcE_step <- function(X.complete, lower, upper, nn, rr, KK, iter, seq.Xm, seq.H, cX, cputime){
  
  
  for(i in 1:nn){
    proc30 = proc.time()[3]
    o = which(rr[i,]==0)
    c = which(rr[i,]!=0)
    ind = seq(from = i, to=KK*nn, by=nn)
    # ==== 1. draw from truncated Univariate Normal Distribution
    if(length(c)==1){
      mp = cond.mean.Precision(mu=seq.Xm[iter,], H=seq.H[[iter]], o=o, c=c, cX=cX[i,])
      
      if(rr[i,c]==-1)
        xs <- tmvtnsim::rtnorm(mean=mp$muc, sd = sqrt(1/mp$Hcc), lower=-Inf,     upper=lower[c], n = KK)
      if(rr[i,c]==1)
        xs <- tmvtnsim::rtnorm(mean=mp$muc, sd = sqrt(1/mp$Hcc), lower=upper[c], upper=Inf,      n = KK)
      X.complete[ind,c] <- xs
      
    }
    proc31 = proc.time()[3]
    cputime["E1"] = cputime["E1"] + proc31 - proc30
    
    # ==== 2. draw from truncated Multivariate Normal Distribution
    if(length(c)>=2){
      LU = LcUc(lower, upper, rr[i,], c); # cat("LU = \n"); print(LU)
      mp = cond.mean.Precision(mu=seq.Xm[iter,], H=seq.H[[iter]], o=o, c=c, cX=cX[i,])
      
      time0 = proc.time()[3]
      # draw samples from truncated Multivariate Normal Distribution
      xs <- tmvtnorm::rtmvnorm.sparseMatrix(n=KK, mean = mp$muc, H=mp$Hcc, 
                                            lower=LU$lower, upper=LU$upper, burn.in.samples=1000, thinning=10)
      X.complete[ind,c] <- xs
      time1 = proc.time()[3]
      cputime["E2"] = cputime["E2"] + time1 - time0
    }
    
  }
  return(list(X.complete=X.complete, cputime=cputime))
}


E_step.exact <- function(cX, rr,lower,upper, mu, H){
  
  # cat("===============begin E_step.exact =========\n")
  pp = ncol(cX)
  nn = nrow(cX)
  
  cX_sub <- cX
  cX_sub[rr!=0] <- 0
  T1_bar <- colSums(cX_sub)
  T2_bar <- crossprod(cX_sub)
  if(any(is.na(T1_bar)) || any(is.na(T2_bar)))
    stop("any(is.na(T1_bar)) || any(is.na(T2_bar))")
  
  for(i in 1:nn){
   
    o = which(rr[i,]==0)
    c = which(rr[i,]!=0)
    
    no <- length(o)
    nc <- length(c)
    x_obs <- cX[i,o]
 
    if(no==pp){
      next
    }
    
    if(nc==pp){
      muc <- mu
      Scc <- solve(H)
    }else{
      mp <-  cond.mean.Precision(mu=mu, H=H, o=o, c=c, cX=cX[i,])
      muc <- mp$muc
      Scc <- solve(mp$Hcc)
    }
  
    LU <- LcUc(lower, upper, rr[i,], c)
    out <-  tmvtnorm::mtmvnorm(mean = as.vector(muc), sigma = Scc, lower =LU$lower, upper=LU$upper)
    tmean <- out$tmean
    tvar <- out$tvar
    
    if(any(is.na(tmean)) || any(is.na(tvar))){
      cat("i=",i,"\n")
      save.image("out.is.na.RData")
      print(out)
      stop("any(is.na(tmean)) || any(is.na(tvar))")
    }  
    
    T1_bar_old = T1_bar
    T2_bar_old = T2_bar
    
   
    T1_bar[c] <- T1_bar[c] + tmean
    T2_bar[o, c] <- T2_bar[o, c, drop = FALSE] + outer(x_obs, tmean)
    T2_bar[c, o] <- t(T2_bar[o, c, drop = FALSE])
    T2_bar[c, c] <- T2_bar[c, c] + tvar + outer(tmean, tmean)
    
    if(any(is.na(T2_bar)) ){
      cat("T2, i=",i,"\n")
      save("T2_bar_old", "T2_bar", "T1_bar_old", "T1_bar",  "out", file="T2.is.na.RData")
      print(T2_bar)
      stop("any(is.na(T2_bar))")
    }  
  }
  
  
  mu_new <- T1_bar / nn
  S_bar <- T2_bar / nn - outer(mu_new, mu_new)
  
  S_bar <- (S_bar + t(S_bar))/2
  T2_bar<- (T2_bar +t(T2_bar))/2
  
  out <- list(mu_new = mu_new, S_bar = S_bar, T1_bar = T1_bar, T2_bar = T2_bar)
  out
}



approx.E_step <- function(cX, rr,lower,upper, mu, H){
  pp = ncol(cX)
  nn = nrow(cX)
 
  cX_sub <- cX
  cX_sub[rr!=0] <- 0
  T1_bar <- colSums(cX_sub)
  T2_bar <- crossprod(cX_sub)
  
  for(i in 1:nn){
    
    o = which(rr[i,]==0)
    c = which(rr[i,]!=0)
  
    no <- length(o)
    nc <- length(c)
    x_obs <- cX[i,o]
   
    if(no==pp){
      next
    }
   
    if(nc==pp){
      muc <- mu
      Scc <- solve(H)
    }else{
      mp <-  cond.mean.Precision(mu=mu, H=H, o=o, c=c, cX=cX[i,])
      muc <- mp$muc
      Scc <- solve(mp$Hcc)
    }
    
    LU <- LcUc(lower, upper, rr[i,], c)
   
    tmean <- vector(mode  = "numeric", length = length(c))
    tvar <- diag(length(c))
    for(m in 1:length(c)){
      out <- tmvtnorm::mtmvnorm(mean = as.vector(muc)[m], sigma = Scc[m, m], lower =LU$lower[m], upper=LU$upper[m])
      tmean[m] <- out$tmean
      tvar[m, m] <- out$tvar
    }
    T1_bar[c] <- T1_bar[c] + tmean
    T2_bar[-c, c] <- T2_bar[-c, c, drop = FALSE] + outer(x_obs, tmean)
    T2_bar[c, -c] <- t(T2_bar[-c, c, drop = FALSE])
    T2_bar[c, c] <- T2_bar[c, c] + tvar + outer(tmean, tmean)
  }
  
  mu_new <- T1_bar / nn
  S_bar <- T2_bar / nn - outer(mu_new, mu_new)
  
  S_bar <- (S_bar + t(S_bar))/2
  T2_bar<- (T2_bar +t(T2_bar))/2
  
  out <- list(mu_new = mu_new, S_bar = S_bar, T1_bar=T1_bar, T2_bar=T2_bar)
  out
}

cEM_ada <- function(cX, rr, lower, upper, lambda = NULL, rho = NULL, L = 10,
                    duplicated = TRUE,
                    E.step = c("exact", "approx", "mcem"), initial=NULL, 
                    penalty = c("L1", "SCAD", "MCP"), ## added for choosing penalty
                    ada_gamma = 3.7, ## added for tuning parameter in SCAD & MCP
                    KK= 200, thr = 1e-03, iter.max=50)
{
  
  pp = ncol(cX)
  nn = nrow(cX)
  
  # ====================== CPU time =====================
  cputime <-rep(0,4)
  names(cputime)<-c("initial", "E-step", "M-step", "total")
  
  proc1 = proc.time()[3]
  
  # ====================== intermediate parameter =======
  seq.H <- seq.S <-seq.Sig<- list()
  seq.Xm<- matrix(0, nrow=iter.max, ncol=pp)
  loglik <-bic<-aic<-loss_mean<-fit_loglik<-loglikpen<- rep(Inf, iter.max)
  # X.complete is used to store the complete samples in E-step
  if(E.step == "mcem"){
    X.complete = cX 
    if(KK >= 2){
      for(i in 2:KK){
        X.complete = rbind(X.complete, cX)
      }
    }
    nK = nrow(X.complete)
  }
  
  # ====================== initial values ===============
  
  seq.Xm[1,]<-colMeans(cX)
  seq.H[[1]]<-solve(cov(cX)+diag(1e-5, nrow = pp))
  seq.S[[1]]  <- cov(cX) #
  seq.Sig[[1]]<-cov(cX) #
 
  df=((sum(seq.Sig[[1]]!=0)-nrow(seq.Sig[[1]]))/2+nrow( seq.Sig[[1]]))
  loglik[1]<- 1
  loglikpen[1]<-loglik[1]+sum(lambda*abs(seq.Sig[[1]]))
  bic[1] <- -2*loglik[1] + log(nn)*df
  aic[1] <- -2*loglik[1] + 2*df
  loss_mean[1]<-0
  # ====================== em iteration ===============
  proc2 = proc.time()[3]
  
  cputime["initial"] = proc2 - proc1
  
  for(iter in 1:(iter.max-1)){
    # cat("iter =", iter, "\n")
    proc3 = proc.time()[3]

    #========== E-step: compute mean and covariance of truncated multivariate normal distribution by R package tmvtnorm::mtmvnorm 
    if(E.step=="exact"){
      result <- E_step.exact(cX, rr,lower,upper, mu=seq.Xm[iter,], H=seq.H[[iter]])
      Xmu <- result$mu_new
      S.cov <- result$S_bar
      
      
      if(det(S.cov)<=0){
        print(S.cov)
        print(eigen(S.cov))
        stop("det(S.cov)<=0.")
      }
      
    }    else if(E.step == "approx"){
      result <- approx.E_step(cX, rr,lower,upper, mu=seq.Xm[iter,], H=seq.H[[iter]])
      
      Xmu <- result$mu_new
      S.cov <- result$S_bar
    }    else if(E.step == "mcem"){
      
      result <- mcE_step(X.complete, lower, upper, nn, rr, KK, iter, seq.Xm, seq.H, cX, cputime)
      cputime    <- result$cputime
      X.complete <-result$X.complete
      
      Xmu <- colMeans(X.complete) 
     
      S.cov <- ((nK-1)/nK)*cov(X.complete)+diag(1e-5, nrow = pp)
      
    }    else{
      stop("E.step must be one char of 'exact', 'approx', and 'mcem'.")
    }
    proc4 = proc.time()[3]
    cputime["E-step"] = cputime["E-step"] + proc4 - proc3
    
    # #========== M-step update parameter
    seq.Xm[iter+1,] = Xmu 
  
    S.cov=S.cov+diag(1e-5, nrow = pp)####
    
    # ---- added for SCAD & MCP ----
    current_sig_abs <- abs(seq.Sig[[iter]])
    
    lambda_array <- array(data = lambda, dim = c(pp,pp,1))

    if(penalty == "SCAD"){
      for(j__ in 1:pp){
        for(k__ in 1:pp){
          if(current_sig_abs[j__,k__] > ada_gamma*lambda){
            lambda_array[j__,k__,1] <- 0
          }
          else if(current_sig_abs[j__,k__] > lambda){
            lambda_array[j__,k__,1] <- (ada_gamma*lambda - current_sig_abs[j__,k__])/(ada_gamma - 1)
          }
        }
      }
    }
    else if(penalty == "MCP"){
      for(j__ in 1:pp){
        for(k__ in 1:pp){
          if(current_sig_abs[j__,k__] >= ada_gamma*lambda){
            lambda_array[j__,k__,1] <- 0
          }
          else{
            lambda_array[j__,k__,1] <- lambda - current_sig_abs[j__,k__]/ada_gamma
          }
        }
      }
    }
    
    fit <- covglasso_ada(S = S.cov, n = nn, lambda = lambda_array,
                         lam = lambda, pen_type = penalty, ga = ada_gamma,
                         rho = rho,
                     )
    # ---- end for added ----
    
    proc7 = proc.time()[3]
    cputime["M-step"] = cputime["M-step"] + proc7 - proc4
    
    seq.Sig[[iter+1]] <- fit$sigma
    seq.S[[iter+1]] = S.cov
    seq.H[[iter+1]] = solve(fit$sigma)
    fit_loglik[[iter+1]]=fit$loglik
    
    loglik[iter+1] <- compute.loglik(n=nn, S=S.cov, wi=seq.H[[iter+1]])
    loglikpen[iter+1]<-loglik[iter+1]+sum(lambda*abs(seq.Sig[[iter+1]]))
    df[iter+1]=((sum(seq.Sig[[iter+1]]!=0)-nrow(fit$sigma))/2+nrow(fit$sigma))
    bic[iter+1] <- -2*loglik[iter+1] + log(nn)*df[iter+1]
    aic[iter+1] <- -2*loglik[iter+1] + 2*df[iter+1]
    
    if(abs((loglikpen[iter+1]-loglikpen[iter])/loglikpen[iter]) < thr){
      break
    }
  }
  
  
  # ==================
  seq = list(#Xm=seq.Xm[1:(iter+1),], 
    #S=seq.S, 
    H=seq.H,
    seq.Sig=seq.Sig
  )
  proc8 = proc.time()[3]
  cputime["total"] = proc8 - proc1
  # ==============输出==========
  if(E.step == "mcem"){
    result = list(Xmu=Xmu, sigma=fit$sigma, H=solve(fit$sigma),
                  fit_loglik=fit_loglik,loglik=loglik,loglikpen=loglikpen,
                  BIC=bic,bic_last=bic[iter+1],AIC=aic,aic_last=aic[iter+1],
                  loglik_last=loglik[iter+1],df=df[iter+1],
                  iter=iter,  cputime=cputime,X.complete=X.complete)
    
  }
  else{
    result = list(Xmu=Xmu, sigma=fit$sigma, H=solve(fit$sigma),
                  fit_loglik=fit_loglik,loglik=loglik,loglikpen=loglikpen,
                  BIC=bic,bic_last=bic[iter+1],df=df[iter+1],
                  AIC=aic,aic_last=aic[iter+1],loglik_last=loglik[iter+1],
                  iter=iter,  cputime=cputime)
    
  }
  
  return(result)
  
}
