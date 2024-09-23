# Generic make method
setGeneric("make_tmb_lists", function(object) standardGeneric("make_tmb_lists"))
# Generic build method
setGeneric("build_tmb", function(object) standardGeneric("build_tmb"))

#' build
#'
#' Create list of tmb data and parameter objects.
#'
#' @param object A data object.
#' @return The fitted tmb model
#' @export
# Method building lists of data and parameters for the model
setMethod("make_tmb_lists", "tmb_list", function(object) {
  frm <- formula(object@formula)

  lm.form <- tryCatch(lme4::lFormula(frm, object@raw_data),
                   error = function(e) e,
                   warning = function(w) w)
  if(length(lm.form$message)>0){
    print(lm.form$message)
    lm.form <- tryCatch(model.matrix(frm, object@raw_data),
                     error = function(e) e,
                     warning = function(w) w)
    model_frame <- model.frame(frm, object@raw_data)
    response_var <- model_frame[[1]]
    lm.form <- list(X = lm.form,
                    fr = response_var)

    object@lm.form <- lm.form
    object@lm.form$model_type <- "lm"
    object@parameters <- list(b = rep(0,ncol(lm.form$X)))
  }else{
    object@lm.form <- lm.form
    object@lm.form$model_type <- "lme"
    object@parameters <- list(b = rep(0,ncol(lm.form$X)),
                              re = rep(0,nrow(lm.form$reTrms$Zt)))
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
  if("reTrms" %in% names(object@lm.form)){
      tmb.data[['y']] <<- stats::model.response(object@lm.form$fr)
      tmb.data[['Zt']] <<- object@lm.form$reTrms$Zt
      tmb.data[['model_type']] <<- object@lm.form$model_type
  }else{
    tmb.data[['y']] <<- object@lm.form$fr
    tmb.data[['Zt']] <<- NA
    tmb.data[['model_type']] <<- object@lm.form$model_type
  }
  tmb.data[['mat']] <<- as.matrix(object@lm.form$X)

  parameters <<- list()
  # parameters[['b']] <<- object@parameters[[1]]
  for(i in names(object@parameters)){
    parameters[[i]] <<- object@parameters[[i]]
  }
  environment(f) <- .GlobalEnv

  if(tmb.data$model_type=="lme"){
    obj <- RTMB::MakeADFun(f,
                           random = c('re'),
                           parameters)
  }else{
    obj <- RTMB::MakeADFun(f,
                           parameters)
  }

  opt <- nlminb(obj$par, obj$fn, obj$gr)

  object@TMB$parameter <- parameters
  object@TMB$tmb.data <- tmb.data
  object@TMB$obj <- obj
  object@TMB$opt <- opt
  object@TMB$sd <- RTMB::sdreport(obj)
  object@TMB$rep <- obj$report()

  return(object)  # Return the updated object with results
})


setMethod("plot", signature(x = "tmb_list", y = "missing"),
          function(x, columns = NULL, nr = 20, marginal = NULL, ...) {
            plot_output <- list()
            # Assuming the original poly() transformation is stored in the S4 object
            v <- x@raw_data[,columns]


            # Extract the terms object
            f <- formula(tmb_model@formula)
            term_obj <- terms(f)

            # Extract the right-hand side (predictor variables) as a character vector
            predictor_terms <- attr(term_obj, "term.labels")
            var_form <- predictor_terms[grep(columns,predictor_terms)]
            var_transform <- eval(parse(text = var_form), x@raw_data)
            # print(var_transform)

            pred_df <- matrix(rep(colMeans(x@TMB$tmb.data$mat), each = 20),
                              nrow = 20,
                              ncol = ncol(x@TMB$tmb.data$mat),
                              byrow = FALSE)

            # Apply poly transformation to new data
            var_seq <- seq(min(v),max(v), length.out=nr)
            new_var <- predict(var_transform, var_seq)

            pred_df[,grep(columns,colnames(x@TMB$tmb.data$mat))] <- new_var
            pred_df <- as.matrix(pred_df)
            pars <- tmb_model@TMB$opt$par

            # print(head(pred_df))
            b <- pars[names(pars)%in%c("b")]
            variances <- diag(pred_df %*% x@TMB$sd$cov.fixed %*% t(pred_df))
            # print(dim(variances))
            # Standard error is the square root of the variance
            standard_error <- sqrt(variances)
            est <- as.matrix(pred_df)%*%b
            upr <- est + 1.96 * standard_error
            lwr <- est - 1.96 * standard_error

            plot_df <- data.frame(est, lwr , upr)

            p <- plot_df %>%
              ggplot2::ggplot(aes (x = var_seq,
                                   y = est), fill = "grey") +
              geom_line() +
              ylab("") +
              xlab(columns) +
              geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1) +
              theme_classic()

            plot_output$p <- p
            plot_output$plot_df <- plot_df
            print(p)
            return(invisible(plot_output))
          })

