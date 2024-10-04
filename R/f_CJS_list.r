f_CJS_list <- function(parms){

  RTMB::getAll(tmb.data,
               parms)

  nll <- 0
  # pred <- rep(0,nrow(p.mat))
  p <- matrix(0,nrow(ch),ncol(ch))
  phi <- matrix(0,nrow(ch),ncol(ch))
  lam <- rep(0,ncol(ch))

  for(i in 1:length(lam.list)){
    lam[i] <- lam.list[[i]]%*%lam.list.b
  }
  if(lam.model_type=="lme"){
    lam <- lam + t(lam.Zt) %*% lam.re
  }


  for(i in 1:length(p.list)){
    for(j in 1:nrow(p.list[[i]])){
      p[i,j+1] <- p.list[[i]][j,]%*%p.list.b
      if(!is.null(dim(p.Zt))){
        p[i,j+1] <- p[i,j+1] + t(p.Zt[,i]) %*% p.re
      }
    }
  }
  p[,ncol(ch)] <- lam
  p <- RTMB::plogis(p)


  for(i in 1:length(phi.list)){
    for(j in 1:nrow(phi.list[[i]])){
      phi[i,j+1] <- phi.list[[i]][j,]%*%phi.list.b
      if(phi.model_type=="lme"){
        phi[i,j+1] <- phi[i,j+1] + t(phi.Zt[,i]) %*% phi.re
      }
    }
  }
  phi[,ncol(ch)] <- lam
  phi <- RTMB::plogis(phi)


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
    nll <- nll - sum(RTMB::dnorm(p.re,0,exp(p.re.sig[p.Lind]),log=TRUE))
    RTMB::REPORT(p.re)
  }
  if(phi.model_type=="lme"){
    nll <- nll - sum(RTMB::dnorm(phi.re,0,exp(phi.re.sig[phi.Lind]),log=TRUE))
    RTMB::REPORT(phi.re)
  }
  if(lam.model_type=="lme"){
    nll <- nll - sum(RTMB::dnorm(lam.re,0,exp(lam.re.sig[lam.Lind]),log=TRUE))
    RTMB::REPORT(lam.re)
  }

  RTMB::REPORT(ch)
  RTMB::REPORT(phi)
  RTMB::REPORT(chi)
  RTMB::REPORT(p)

  return(nll)
}

# tmb_obj@TMB
