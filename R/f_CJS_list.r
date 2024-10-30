f <- function(parms){

  RTMB::getAll(tmb.data,
               parms)

  nll <- 0
  neglls <- list()
  # pred <- rep(0,nrow(p.mat))
  state_dim <- length(MR_settings$states)
  p <- array(Inf,c(state_dim-1,nrow(ch),ncol(ch)))
  phi <- array(Inf,c(state_dim-1,nrow(ch),ncol(ch)))
  eta <- matrix(-Inf,nrow(ch),ncol(ch))
  lam <- matrix(0,state_dim-1,nrow(ch))


  #For each capture, calculate the probability based on gamma and omega
  #Different capture history for each location
  gamma <- array(0,c(nrow(ch),ncol(ch),length(states),length(states)),
                 dimnames = list(id = seq_len(nrow(ch)),
                                 loc = names(ch),
                                 init_state = states,
                                 next_state = states))
  omega <- array(0,c(nrow(ch),ncol(ch),length(states),length(states)),
                 dimnames = list(id = seq_len(nrow(ch)),
                                 loc = names(ch),
                                 init_state = states,
                                 next_state = states))


  for(state in 1:(state_dim-1)){
    for(i in 1:length(lam.list)){
      lam[state,i] <- lam.list[[i]]%*%lam.list.b#[state,]
    }
    if(lam.model_type=="lme"){
      lam[state,] <- lam[state,] + t(lam.Zt) %*% lam.re
    }
  }

  for(state in 1:(state_dim-1)){
    icnt <- 1
    for(i in 1:length(p.list)){
      for(j in 1:nrow(p.list[[i]])){
        p[state,i,j+1] <- p.list[[i]][j,]%*%p.list.b[state,]
        if(!is.null(dim(p.Zt))){
          p[state,i,j+1] <- p[state,i,j+1] + t(p.Zt[,icnt]) %*% p.re
          icnt <- icnt + 1
        }
      }
    }
    p[state,,ncol(ch)] <- lam[state,]
    # print(p[state,,ncol(ch)])
  }
  p <- RTMB::plogis(p)


  for(state in 1:(state_dim-1)){
    icnt <- 1
    for(i in 1:length(phi.list)){
      for(loc in 1:nrow(phi.list[[i]])){
        phi[state,i,loc+1] <- phi.list[[i]][loc,]%*%phi.list.b[state,]
      }
      if(phi.model_type=="lme"){
        phi[state,i,loc+1] <- phi[state,i,loc+1] + t(phi.Zt[,icnt]) %*% phi.re
        icnt <- icnt + 1
      }
    }
    phi[state,,ncol(ch)] <- lam[state,]
    # print(phi[state,,ncol(ch)])
  }
  phi <- RTMB::plogis(phi)


  if(state_dim>2 & length(eta.list)>1){
    icnt <- 1
    for(i in 1:length(eta.list)){
      for(loc in 1:nrow(eta.list[[i]])){
        eta[i,loc+1] <- eta.list[[i]][loc,]%*%eta.list.b
        if(eta.model_type=="lme"){
          eta[i,loc+1] <- eta[i,loc+1] + t(eta.Zt[,icnt]) %*% eta.re
          icnt <- icnt + 1
        }
      }
    }
    eta <- RTMB::plogis(eta)
  }else{
   # eta <- 0#RTMB::plogis(eta)
  }

  if(model_type == "CJS"){
    chi <- rep(0, nrow(ch))
    for(i in 1:nrow(ch)){
      chi[i] <- 1
      if(U[i]<ncol(ch)){
        for(j in (ncol(ch)-1):U[i]){
          chi[i] <- (1 - phi[1,i,j+1]) + phi[1,i,j+1] * (1 - p[1,i,j+1]) * chi[i]
        }
      }
    }

    for(i in 1:nrow(ch)){
      if(U[i]>1){
        for(j in 2:U[i]){
          nll <- nll - RTMB::dbinom(1, 1, phi[1,i,j],log=TRUE) * n[i]
          nll <- nll - RTMB::dbinom(ch[i,j], 1, p[1,i,j],log=TRUE)  * n[i]
        }
      }
      nll <- nll - RTMB::dbinom(1,1,chi[i],log=TRUE) * n[i]
    }
  }else{
    #Collect the negative log likelihoods while running the forward algorithm
    ch[ch==0] <- max(ch) + 1
    max_s <- max(ch)
    for (i in seq_len(nrow(ch))) {

      #Set initial value of fish based on pool where it was tagged and released
      delta <- rep(0,max_s)
      delta[init_state[i]] <- 1 #delta[1] <- 1 means subYr, delta[2] <- 1 means yr

      #t is a dam location not time as in the case of Labuzzeta
      for(loc in 2:ncol(ch)){
        for(state in 1:(max_s-1)){
          gamma[i,loc,state,state] <- RTMB::qlogis(phi[state,i,loc])
          if(state_dim>2){
            if(state==1){
              gamma[i,loc,state,state+1] <- RTMB::qlogis(eta[i,loc])
            }else{
              gamma[i,loc,state,state+1] <- -Inf
            }
          }
          gamma[i,loc,state,max_s] <- 1
          gamma[i,loc,state,state:max_s] <- exp(gamma[i,loc,state,state:max_s])/sum(exp(gamma[i,loc,state,state:max_s]))
          omega[i,loc,state,state] <- p[state,i,loc]
          omega[i,loc,state,max_s] <- 1 - p[state,i,loc]

          }
          # print(paste(i, loc))
        # print(gamma[loc,,])
        gamma[i,loc,max_s,max_s] <- 1
        omega[i,loc,max_s,max_s] <- 1
      }

      Ui <- init_loc[i] #Start at first location
      prods <- list(diag(1, max_s))
      for (loc in (Ui + 1):ncol(ch)) {
        prods[[loc - Ui + 1]] <- prods[[loc - Ui]] %*% gamma[i,loc,,] %*% RTMB::diag(omega[i,loc,, ch[i,loc]])
      }
      #Keep the likelihood a "true" likelihood for simulation purposes
      nll <- nll - RTMB::dbinom(1,1,sum(t(delta) %*% prods[[length(prods)]]), log = TRUE) * n[i]
    }
  }

  if(p.model_type=="lme"){
    nll <- nll - sum(RTMB::dnorm(p.re, 0, exp(p.re.sig[p.Lind]), log=TRUE))
    RTMB::REPORT(p.re)
  }
  if(phi.model_type=="lme"){
    nll <- nll - sum(RTMB::dnorm(phi.re, 0, exp(phi.re.sig[phi.Lind]), log=TRUE))
    RTMB::REPORT(phi.re)
  }
  if(lam.model_type=="lme"){
    nll <- nll - sum(RTMB::dnorm(lam.re, 0, exp(lam.re.sig[lam.Lind]), log=TRUE))
    RTMB::REPORT(lam.re)
  }
  if(eta.model_type=="lme"){
    nll <- nll - sum(RTMB::dnorm(eta.re, 0, exp(eta.re.sig[eta.Lind]), log=TRUE))
    RTMB::REPORT(eta.re)
  }


  RTMB::REPORT(ch)
  RTMB::REPORT(phi)
  # RTMB::REPORT(chi)
  RTMB::REPORT(p)
  RTMB::REPORT(gamma)
  RTMB::REPORT(omega)
  RTMB::REPORT(lam)
  if(length(eta.list)>1){
    RTMB::REPORT(eta)
  }

  # y_i <- as.integer(factor(raw$Year, levels = sort(unique(raw$Year))))
  # loc_i <- as.integer(factor(raw$loc, levels = locs))
  # rs_i <- as.integer(factor(raw$ReleaseSite))
  #
  # df <- array(0, c(length(unique(y_i)),length(unique(loc_i))+1, length(unique(rs_i))))
  # print(dim(df))
  # # print(list(y_i = levels(y_i), loc_i = c(levels(loc_i),'lam'), rs_i = levels(rs_i)))
  # # dimnames(df) <- list(y_i = sort(unique(raw$Year)), loc_i = locs, rs_i = levels(factor(raw$ReleaseSite))
  #
  # # names(df) <- locs
  # for(i in 1:nrow(raw)){
  #   df[y_i[i],loc_i[i]-1,rs_i[i]] <- df[y_i[i],loc_i[i]-1,rs_i[i]] + prod(phi[i,1:loc_i[i]]) * p[i,loc_i[i]] * n[i]
  # }
  # for(i in 1:nrow(raw)){
  #   df[y_i[i],max(loc_i),rs_i[i]] <- df[y_i[i],max(loc_i),rs_i[i]] + prod(phi[i,1:(max(loc_i)+1)]) * p[i,(max(loc_i)+1)] * n[i]
  # }
  # sm <- aggregate(df[,locs],by=list(year = raw$Year),sum)[,locs[2:length(locs)]]
  # RTMB::REPORT(df)
  # RTMB::REPORT(y_i)
  # RTMB::REPORT(loc_i)
  # RTMB::ADREPORT(df)
  # RTMB::ADREPORT(do.call("C",sm))

  return(nll)
}

# tmb_obj@TMB
