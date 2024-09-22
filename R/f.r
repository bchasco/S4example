f <- function(parms){

  RTMB::getAll(tmb.data,
               parms)

  pred <- rep(0,nrow(mat))
  for(i in seq_along(pred)){
    pred[i] <- sum(mat[i,] * b)
  }
  # f_b <- plogis(b)
  nll <- sum((y-pred)^2)
  RTMB::REPORT(mat)
  return(nll)
}

# tmb_obj@TMB
