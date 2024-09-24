# Generic make method
setGeneric("make_tmb_lists", function(object) standardGeneric("make_tmb_lists"))
# Generic build method
setGeneric("build_tmb", function(object) standardGeneric("build_tmb"))
# Define the generic function first
setGeneric("marginal_plot", function(x, y = NULL, ...) {
  standardGeneric("marginal_plot")
})

#' build
#'
#' Create list of tmb data and parameter objects.
#'
#' @param object A data object.
#' @return The fitted tmb model
#' @export
# Method building lists of data and parameters for the model
setMethod("make_tmb_lists", "tmb_list", function(object) {
  p.frm <- formula(object@p.frm)
  phi.frm <- formula(object@phi.frm)
  lam.frm <- formula(object@lam.frm)

  #Detection probability
  p.lm.form <- tryCatch(lme4::lFormula(p.frm, object@raw_data),
                      error = function(e) e,
                      warning = function(w) w)
  if(length(p.lm.form$message)>0){
    print(p.lm.form$message)
    p.lm.form <- tryCatch(model.matrix(p.frm, object@raw_data),
                        error = function(e) e,
                        warning = function(w) w)
    p.df <- model.frame(p.frm, object@raw_data)
    p.res <- p.df[[1]]
    p.lm.form <- list(X = p.lm.form,
                    fr = p.res)

    object@p.lm.form <- p.lm.form
    object@p.lm.form$p.model_type <- "lm"
    object@parameters <- append(object@parameters,
                                list(p.b = rep(0,ncol(p.lm.form$X))))
  }else{
    object@p.lm.form <- p.lm.form
    object@p.lm.form$p.model_type <- "lme"
    object@parameters <- append(object@parameters,
                                list(p.b = rep(0,ncol(p.lm.form$X)),
                                     p.re = rep(0,nrow(p.lm.form$reTrms$Zt)),
                                     p.re.sig = length(p.lm.form$reTrms$cnms)))
  }

  #survival forms
  phi.lm.form <- tryCatch(lme4::lFormula(phi.frm, object@raw_data),
                   error = function(e) e,
                   warning = function(w) w)
  if(length(phi.lm.form$message)>0){
    print(phi.lm.form$message)
    phi.lm.form <- tryCatch(model.matrix(phi.frm, object@raw_data),
                     error = function(e) e,
                     warning = function(w) w)
    phi.df <- model.frame(phi.frm, object@raw_data)
    phi_res <- phi.df[[1]]
    phi.lm.form <- list(X = phi.lm.form,
                    fr = phi_res)

    object@phi.lm.form <- phi.lm.form
    object@phi.lm.form$phi.model_type <- "lm"
    object@parameters <- append(object@parameters,
                                                     list(phi.b = rep(0,ncol(phi.lm.form$X))))
  }else{
    object@phi.lm.form <- phi.lm.form
    object@phi.lm.form$phi.model_type <- "lme"
    object@parameters <- append(object@parameters,
                                                     list(phi.b = rep(0,ncol(phi.lm.form$X)),
                                                          phi.re = rep(0,nrow(phi.lm.form$reTrms$Zt)),
                                                          phi.re.sig = length(phi.lm.form$reTrms$cnms)))
  }

  #lam forms
  lam.lm.form <- tryCatch(lme4::lFormula(lam.frm, object@raw_data),
                          error = function(e) e,
                          warning = function(w) w)
  if(length(lam.lm.form$message)>0){
    print(lam.lm.form$message)
    lam.lm.form <- tryCatch(model.matrix(lam.frm, object@raw_data),
                            error = function(e) e,
                            warning = function(w) w)
    lam.df <- model.frame(lam.frm, object@raw_data)
    lam_res <- lam.df[[1]]
    lam.lm.form <- list(X = lam.lm.form,
                        fr = lam_res)

    object@lam.lm.form <- lam.lm.form
    object@lam.lm.form$lam.model_type <- "lm"
    object@parameters <- append(object@parameters,
                                list(lam.b = rep(0,ncol(lam.lm.form$X))))
  }else{
    object@lam.lm.form <- lam.lm.form
    object@lam.lm.form$lam.model_type <- "lme"
    object@parameters <- append(object@parameters,
                                list(lam.b = rep(0,ncol(lam.lm.form$X)),
                                     lam.re = rep(0,nrow(lam.lm.form$reTrms$Zt)),
                                     lam.re.sig = length(lam.lm.form$reTrms$cnms)))
  }

  return(object)
})

#' build
#'
#' Builds the tmb objective
#'
#' @param object A data object.
#' @return The tmb objective function
#' @export
# Method building and runnning the tmb model
setMethod("build_tmb", "tmb_list", function(object) {
  tmb.data <<- list() #create in .GlobalEnv
  tmb.data[['n']] <<- object@raw_data$n
  if("reTrms" %in% names(object@phi.lm.form)){
      tmb.data[['phi.Zt']] <<- as.matrix(object@phi.lm.form$reTrms$Zt)
      tmb.data[['phi.model_type']] <<- object@phi.lm.form$phi.model_type
  }else{
    tmb.data[['phi.Zt']] <<- NA
    tmb.data[['phi.model_type']] <<- object@phi.lm.form$phi.model_type
  }
  if("reTrms" %in% names(object@p.lm.form)){
    tmb.data[['p.Zt']] <<- as.matrix(object@p.lm.form$reTrms$Zt)
    tmb.data[['p.model_type']] <<- object@p.lm.form$p.model_type
  }else{
    tmb.data[['p.Zt']] <<- NA
    tmb.data[['p.model_type']] <<- object@p.lm.form$p.model_type
  }
  if("reTrms" %in% names(object@lam.lm.form)){
    tmb.data[['lam.Zt']] <<- as.matrix(object@lam.lm.form$reTrms$Zt)
    tmb.data[['lam.model_type']] <<- object@lam.lm.form$lam.model_type
  }else{
    tmb.data[['lam.Zt']] <<- NA
    tmb.data[['lam.model_type']] <<- object@lam.lm.form$lam.model_type
  }


  tmb.data[['p.mat']] <<- as.matrix(object@p.lm.form$X)
  tmb.data[['phi.mat']] <<- as.matrix(object@phi.lm.form$X)
  tmb.data[['lam.mat']] <<- as.matrix(object@lam.lm.form$X)

  tmb.data[['ch']] <<- as.matrix(object@raw_data[,c('Tagged', 'as.Smolt', 'As.Adult.ballard')])
  tmb.data[["U"]] <<- apply(tmb.data[['ch']],1,function(x){max(which(x>0,arr.ind = TRUE))})

  parameters <<- list()
  # parameters[['b']] <<- object@parameters[[1]]
  for(i in names(object@parameters)){
    parameters[[i]] <<- object@parameters[[i]]
  }
  environment(f_CJS) <- .GlobalEnv

  random <<- c()
  if(tmb.data$p.model_type=="lme"){
    random <<- c(random,"p.re")
  }
  if(tmb.data$phi.model_type=="lme"){
    random <<- c(random,"phi.re")
  }

  obj <- RTMB::MakeADFun(f_CJS,
                         random = random,
                         parameters)

  opt <- nlminb(obj$par, obj$fn, obj$gr)

  object@TMB$parameter <- parameters
  object@TMB$tmb.data <- tmb.data
  object@TMB$obj <- obj
  object@TMB$opt <- opt
  object@TMB$sd <- RTMB::sdreport(obj)
  object@TMB$rep <- obj$report()

  return(object)  # Return the updated object with results
})


setMethod("marginal_plot", signature(x = "tmb_list", y = "missing"),
          function(x, columns = NULL, nr = 20, proc = "p", marginal = NULL, ...) {
            plot_output <- list()
            # Assuming the original poly() transformation is stored in the S4 object
            v <- x@raw_data[,columns]


            # Accessing the slot dynamically
            slot_name <- paste0(proc, ".frm")  # Create the dynamic slot name
            # Extract the terms object
            f <- formula(slot(x, slot_name))
            term_obj <- terms(f)

            # Extract the right-hand side (predictor variables) as a character vector
            predictor_terms <- attr(term_obj, "term.labels")
            var_form <- predictor_terms[grep(columns,predictor_terms)]
            var_transform <- eval(parse(text = var_form), x@raw_data)

            slot_name <- paste0(proc, ".lm.form")  # Create the dynamic slot name
            slot_value <- slot(x, slot_name)  # Access the slot dynamically
            pred_df <- matrix(rep(colMeans(slot_value$X), each = 20),
                              nrow = 20,
                              ncol = ncol(slot_value$X),
                              byrow = FALSE)

            # # Apply poly transformation to new data
            var_seq <- seq(min(v),max(v), length.out=nr)
            # if(grep("poly",))
            new_var <- predict(var_transform, var_seq)
            #
            pred_df[,grep(columns,colnames(slot_value$X))] <- new_var
            pred_df <- as.matrix(pred_df)

            par_idx <- grep(paste0(proc,".b"),names(x@TMB$opt$par))
            pars <- x@TMB$opt$par[par_idx]
            variances <- diag(pred_df %*% x@TMB$sd$cov.fixed[par_idx,par_idx] %*% t(pred_df))

            # # Standard error is the square root of the variance
            standard_error <- sqrt(variances)


            est <- as.matrix(pred_df)%*%pars
            upr <- est + 1.96 * standard_error
            lwr <- est - 1.96 * standard_error

            plot_df <- data.frame(est, lwr , upr)

            p <- plot_df %>%
              ggplot2::ggplot(aes (x = var_seq,
                                   y = plogis(est)), fill = "grey") +
              geom_line() +
              ylab("") +
              xlab(columns) +
              geom_ribbon(aes(ymin = plogis(lwr), ymax = plogis(upr)), alpha = 0.1) +
              theme_classic()

            plot_output$p <- p
            plot_output$plot_df <- plot_df
            print(p)

            return(invisible(plot_output))
          })

