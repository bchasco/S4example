f_CJS <- function(parms){

  RTMB::getAll(tmb.data,
               parms)

  nll <- 0
  pred <- rep(0,nrow(p.mat))
  p <- matrix(0,nrow(ch),ncol(ch))
  phi <- matrix(0,nrow(ch),ncol(ch))
  lam <- rep(0,ncol(ch))

  lam <- lam.mat%*%lam.b
  if(lam.model_type=="lme"){
    lam <- lam + t(lam.Zt) %*% lam.re
  }

  p[,2] <- p.mat%*%p.b
  if(p.model_type=="lme"){
    p[,2] <- p[,2] + t(p.Zt) %*% p.re
  }
  p[,2] <- RTMB::plogis(p[,2])
  p[,3] <- RTMB::plogis(lam)

  phi[,2] <- phi.mat%*%phi.b
  if(phi.model_type=="lme"){
    phi[,2] <- phi[,2] + t(phi.Zt) %*% phi.re
  }
  phi[,2] <- RTMB::plogis(phi[,2])
  phi[,3] <- RTMB::plogis(lam)


  chi <- rep(0, nrow(ch))
  for(i in 1:nrow(ch)){
    chi[i] <- 1
    if(U[i]<ncol(ch)){
      for(j in (ncol(ch)-1):U[i]){
        chi[i] <- (1 - phi[i,j+1]) + phi[i,j+1] * (1 - p[i,j+1]) * chi[i]
      }
    }
  }

  for(i in 1:nrow(ch)){
    if(U[i]>1){
      for(j in 2:U[i]){
        nll <- nll - RTMB::dbinom(1,1,phi[i,j],log=TRUE) * n[i]
        nll <- nll - RTMB::dbinom(ch[i,j],1,p[i,j],log=TRUE)  * n[i]
      }
    }
    nll <- nll - RTMB::dbinom(1,1,chi[i],log=TRUE) * n[i]
  }

  if(p.model_type=="lme"){
    nll <- nll - sum(RTMB::dnorm(p.re,0,exp(p.re.sig),log=TRUE))
    RTMB::REPORT(p.re)
  }
  if(phi.model_type=="lme"){
    nll <- nll - sum(RTMB::dnorm(phi.re,0,exp(phi.re.sig),log=TRUE))
    RTMB::REPORT(phi.re)
  }
  if(lam.model_type=="lme"){
    nll <- nll - sum(RTMB::dnorm(lam.re,0,exp(lam.re.sig),log=TRUE))
    RTMB::REPORT(lam.re)
  }

  RTMB::REPORT(p.mat)
  RTMB::REPORT(phi.mat)
  RTMB::REPORT(ch)
  RTMB::REPORT(phi)
  RTMB::REPORT(chi)
  RTMB::REPORT(p)

  # if(model_type=="lme"){
  #   RTMB::REPORT(re)
  # }
  return(nll)
}

# tmb_obj@TMB
