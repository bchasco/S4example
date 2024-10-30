f <- function(parms){

  RTMB::getAll(tmb.data,
               parms)

  nll <- 0
  pred <- rep(0,nrow(mat))
  for(i in seq_along(pred)){
    pred[i] <- sum(mat[i,] * b)
    if(model_type=="lme"){
      pred[i] <- pred[i] + sum(re * Zt[,i])
    }
  }
  if(model_type=="lme"){
    nll <- nll - sum(log(RTMB::dnorm(re,0,10)))
  }

  nll <- nll + sum((y-pred)^2)

  RTMB::REPORT(mat)
  if(model_type=="lme"){
    RTMB::REPORT(re)
  }
  return(nll)
}

# tmb_obj@TMB
