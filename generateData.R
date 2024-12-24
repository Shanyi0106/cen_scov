generate.true.graph.Para.Data<-function(nn, pp, Sigma, prob_cens,
                                        lower=lower, 
                                        upper=rep(Inf,pp)){
  
  # =======================================================
  # 1. generate true parameter values
  # =======================================================
  # pp  # the nubmer of variables in Gaussian distribution
  
  #============== 1.2  parameters 
  mu = rep(0, pp) # mean
  
  tau<-qnorm(prob_cens)
  lower=sqrt(diag(Sigma))*tau
  HH =  solve(Sigma)
  
  # =======================================================
  # 2. draw samples from censored Gaussian
  # =======================================================
  # =========== 2.1 draw samples from multivariate normal distribution
  
  XX = MASS::mvrnorm(n=nn, mu=mu, Sigma = Sigma)
  colnames(XX) = paste("V", 1:pp, sep="")
  
  # =========== 2.2 draw samples from censored multivariate normal distribution
  cX = XX 
  
  rr = matrix(0, nrow=nn, ncol=pp)
  for(i in 1:pp){
    # censored below
    index = which(XX[,i]<lower[i])
   
    cX[index, i] = lower[i]
    rr[index, i] = -1
    
    # censored above
    index = which(XX[,i]>upper[i])
    cX[index, i] = upper[i]
   
    rr[index, i] = 1
  }
  
  # compute stat for cX 
  sum(rr==1)/(nn*pp)
  sum(rr==0)/(nn*pp)
  sum(rr==-1)/(nn*pp)
  
  co = rowSums(abs(rr))
  length(which(co>0))/nn
  
  truePara = list(mu=mu,Sigma=Sigma, H=HH, 
                  lower=lower, upper=upper)
  list(cX=cX, rr=rr, XX=XX,
       truePara=truePara)
}