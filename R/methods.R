# Generic make method
setGeneric("make_tmb_lists", function(object) standardGeneric("make_tmb_lists"))
# Generic build method
setGeneric("build_tmb", function(object) standardGeneric("build_tmb"))
# Define the generic function first

setGeneric("AIC", function(object) standardGeneric("AIC"))
# Define the generic function first

setGeneric("plot", function(obj, y = NULL, ...) {
  standardGeneric("plot")
})
#AIC calculation
setGeneric("AIC", function(object) standardGeneric("AIC"))


#' build
#'
#' Create list of tmb data and parameter objects.
#'
#' @param object A data object.
#' @return The fitted tmb model
#' @export
# Method building lists of data and parameters for the model
setMethod("make_tmb_lists", "tmb_list", function(object) {
  for(i in c('p','phi','lam')){

    if(i!='lam'){
      state_loc <- 'as.Smolt'
    }else{
      state_loc <- 'As.Adult.ballard'
    }

    i.lm.form <- paste0(i,'.lm.form')

    assign(paste0(i,'.loc'),
           object@raw_data %>%
            mutate(id = row_number()) %>%
            arrange(id) %>%
            pivot_longer(cols = c(Tagged, as.Smolt, As.Adult.ballard),
                         names_to = 'loc',
                         values_to = 'state') %>%
            filter(loc %in% state_loc))

    assign(i.lm.form,
           tryCatch(lme4::lFormula(formula(object@MR_settings$frms[i][[1]]), get(paste0(i,'.loc'))),
                    error = function(e) e,
                    warning = function(w) w)
           )

    if(length(get(i.lm.form)$message)>0){
      #No random-effects
      assign(i.lm.form,
             tryCatch(model.matrix(formula(object@MR_settings$frms[i][[1]]), get(paste0(i,'.loc'))),
                      error = function(e) e,
                      warning = function(w) w))
      assign(paste0(i,'.mat'),
                    as.data.frame(get(i.lm.form)) %>%
                      mutate(id = get(paste0(i,'.loc'))$id) %>%
                      group_by(id) %>%
                      summarise(y_values = list(pick(everything()))))

      assign(i.lm.form, list(X = get(i.lm.form)))
      slot(object,i.lm.form)[[paste0(i,'.data')]] <- get(paste0(i,'.loc'))
      slot(object,i.lm.form)[[paste0(i,'.mat')]] <- get(paste0(i,'.mat'))[[2]]
      slot(object,i.lm.form)[[i.lm.form]] <- get(i.lm.form)
      slot(object,i.lm.form)[[paste0(i,".model_type")]] <- "lm"

    }else{
      #Mixed-effect model
      slot(object,i.lm.form)[[i.lm.form]] <- get(i.lm.form)
      slot(object,i.lm.form)[[paste0(i,".model_type")]] <- "lme"
      slot(object,i.lm.form)[[paste0(i,'.data')]] <- get(paste0(i,'.loc'))
      assign(paste0(i,'.mat'),
             as.data.frame(slot(object,i.lm.form)[[i.lm.form]]$X) %>%
               mutate(id = get(paste0(i,'.loc'))$id) %>%
               group_by(id) %>%
               summarise(y_values = list(pick(everything()))))
      slot(object,i.lm.form)[[paste0(i,'.mat')]] <- get(paste0(i,'.mat'))[[2]]
    }

    if(slot(object, i.lm.form)[[paste0(i,".model_type")]] == "lm"){
      slot(object,'parameters')[[paste0(i,'.list.b')]] <- rep(0,ncol(slot(object,i.lm.form)[[paste0(i,'.mat')]][[1]]))
    }else{
      slot(object,'parameters')[[paste0(i,'.list.b')]] <- rep(0,ncol(slot(object,i.lm.form)[[paste0(i,'.mat')]][[1]]))
      slot(object,'parameters')[[paste0(i,'.re')]] <- rep(0,nrow(slot(object,i.lm.form)[[i.lm.form]]$reTrms$Zt))
      slot(object,'parameters')[[paste0(i,'.re.sig')]] <- rep(0,length(slot(object,i.lm.form)[[i.lm.form]]$reTrms$cnms))
    }

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
  if("reTrms" %in% names(object@phi.lm.form$phi.lm.form)){
      tmb.data[['phi.Zt']] <<- as.matrix(object@phi.lm.form$phi.lm.form$reTrms$Zt)
  }else{
    tmb.data[['phi.Zt']] <<- NA
  }
  tmb.data[['phi.model_type']] <<- object@phi.lm.form$phi.model_type

  if("reTrms" %in% names(object@p.lm.form$p.lm.form)){
    tmb.data[['p.Zt']] <<- as.matrix(object@p.lm.form$p.lm.form$reTrms$Zt)
  }else{
    tmb.data[['p.Zt']] <<- NA
  }
  tmb.data[['p.model_type']] <<- object@p.lm.form$p.model_type

  if("reTrms" %in% names(object@lam.lm.form$lam.lm.form)){
    tmb.data[['lam.Zt']] <<- as.matrix(object@lam.lm.form$lam.lm.form$reTrms$Zt)
  }else{
    tmb.data[['lam.Zt']] <<- NA
  }
  tmb.data[['lam.model_type']] <<- object@lam.lm.form$lam.model_type



  tmb.data[['ch']] <<- as.matrix(object@raw_data[,c('Tagged', 'as.Smolt', 'As.Adult.ballard')])
  tmb.data[["U"]] <<- apply(tmb.data[['ch']],1,function(x){max(which(x>0,arr.ind = TRUE))})

  tmb.data[['phi.list']] <<- list()
  for(i in 1:length(object@phi.lm.form$phi.mat)){
    tmb.data[['phi.list']][[i]] <<- as.matrix(object@phi.lm.form$phi.mat[[i]])
  }

  tmb.data[['lam.list']] <<- list()
  for(i in 1:length(object@lam.lm.form$lam.mat)){
    tmb.data[['lam.list']][[i]] <<- as.matrix(object@lam.lm.form$lam.mat[[i]])
  }

  tmb.data[['p.list']] <<- list()
  for(i in 1:length(object@p.lm.form$p.mat)){
    tmb.data[['p.list']][[i]] <<- as.matrix(object@p.lm.form$p.mat[[i]])
  }

  parameters <<- list()
  # parameters[['b']] <<- object@parameters[[1]]
  for(i in names(object@parameters)){
    parameters[[i]] <<- object@parameters[[i]]
  }
  environment(f_CJS_list) <- .GlobalEnv

  random <<- c()
  if(tmb.data$p.model_type=="lme"){
    random <<- c(random,"p.re")
  }
  if(tmb.data$phi.model_type=="lme"){
    random <<- c(random,"phi.re")
  }
  if(tmb.data$lam.model_type=="lme"){
    random <<- c(random,"lam.re")
  }

  obj <- RTMB::MakeADFun(f_CJS_list,
                         random = random,
                         parameters)

  opt <- nlminb(obj$par, obj$fn, obj$gr)

  object@TMB$parameter <- parameters
  object@TMB$tmb.data <- tmb.data
  object@TMB$obj <- obj
  object@TMB$opt <- opt
  if(object@MR_settings$make_sd){
    object@TMB$sd <- RTMB::sdreport(obj)
  }
  object@TMB$rep <- obj$report()

  return(object)  # Return the updated object with results
})


setMethod("plot", signature(obj = "tmb_list", y = "missing"),
          function(obj, var = NULL, factor = NULL, nr = 20, process = "p", type = NULL, ...) {
            plot_output <- list()

            # Helper function to create prediction data frame
            create_pred_df <- function(slot_value, var_seq, var, new_var, nr) {
              pred_df <- matrix(rep(colMeans(slot_value$X), each = nr),
                                nrow = nr, ncol = ncol(slot_value$X), byrow = FALSE)
              pred_df[, grep(var, colnames(slot_value$X))] <- new_var
              return(as.matrix(pred_df))
            }

            # Helper function to plot the result
            plot_response <- function(var_seq, est, lwr, upr, var_label, y_label) {
              plot_df <- data.frame(est, lwr, upr)
              p <- plot_df %>%
                ggplot2::ggplot(aes(x = var_seq, y = plogis(est)), fill = "grey") +
                geom_line() +
                xlab(var_label) +
                ylim(0, 1) +
                geom_ribbon(aes(ymin = plogis(lwr), ymax = plogis(upr)), alpha = 0.1) +
                theme_classic() +
                ylab(y_label)
              return(list(p = p, plot_df = plot_df))
            }

            par_names <- names(obj@TMB$opt$par)

            if (type == "response") {
              v <- obj@raw_data[, var]
              slot_name <- paste0("MR_settings")
              f <- slot(obj, slot_name)
              term_obj <- terms(formula(f$frms[process][[1]]))
              predictor_terms <- attr(term_obj, "term.labels")
              var_form <- predictor_terms[grep(var, predictor_terms)]
              var_transform <- eval(parse(text = var_form), obj@raw_data)

              slot_value <- slot(obj, paste0(process, ".lm.form"))
              var_seq <- seq(min(v), max(v), length.out = nr)
              new_var <- predict(var_transform, var_seq)

              pred_df <- create_pred_df(slot_value, var_seq, var, new_var, nr)
              par_idx <- grep(paste0(process, ".list.b"), names(obj@TMB$opt$par))
              pars <- obj@TMB$opt$par[par_idx]
              variances <- diag(pred_df %*% obj@TMB$sd$cov.fixed[par_idx, par_idx] %*% t(pred_df))
              standard_error <- sqrt(variances)
              est <- pred_df %*% pars
              upr <- est + 1.96 * standard_error
              lwr <- est - 1.96 * standard_error

              y_label <- ifelse(var == "phi", "Survival probability", "Detection probability")
              result <- plot_response(var_seq, est, lwr, upr, var, y_label)
              plot_output$p <- result$p
              plot_output$plot_df <- result$plot_df
            }

            if (type == "re") {
              int <- obj@TMB$opt$par[par_names%in%paste0(process, '.list.b')][1]
              se_int <- diag(obj@TMB$sd$cov.fixed)[par_names %in% paste0(process, '.list.b')][1]

              re_list <- obj@TMB$obj$env$parList()[names(obj@TMB$obj$env$parList()) %in% obj@TMB$obj$env$.random]
              re_names <- rep(names(re_list), sapply(re_list, length))
              est <- re_list[[paste0(process, '.re')]] + int
              se <- obj@TMB$sd$diag.cov.random[re_names %in% paste0(process, '.re')]
              upr <- est + 1.96 * (se + se_int)
              lwr <- est - 1.96 * (se + se_int)

              plot_df <- data.frame(est, lwr, upr)
              p <- plot_df %>%
                ggplot2::ggplot(aes(x = 1:length(est), y = plogis(est)), fill = "grey") +
                geom_point() +
                ylab(process) +
                ylim(0, 1) +
                geom_errorbar(aes(ymin = plogis(lwr), ymax = plogis(upr)), alpha = 0.5) +
                theme_classic()

              if (process == "phi") {
                p <- p + ylab("Survival probability")
              }
              if (process == "p") {
                p <- p + ylab("Detection probability")
              }
              plot_output$p <- p
              plot_output$plot_df <- plot_df
            }

            print(plot_output$p)
            return(invisible(plot_output))
          })


#' AIC
#'
#' Calculates AIC for tmb_list
#'
#' @param object A tmb list.
#' @return The AIC
#' @export
# Method building and runnning the tmb model
setMethod("AIC", "tmb_list", function(object) {
  n <- sum(object@TMB$opt$n)
  K <- length(object@TMB$opt$par)
  nll <- object@TMB$opt$objective
  AICc <- 2 * nll + 2 * K +(2 * K * ( K + 1) )/(n - K  - 1)
  return(AICc)  # Return the updated object with results
})
