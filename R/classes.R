# Define the general Model class
setClass(
  "Model",
  slots = list(
    data = "ANY",
    model_type = "character",
    fit_result = "ANY"
  )
)

#' tmb_data Class
#'
#' This class holds the data.
#'
#' @slot raw_data A list of data.
#' @slot TMB
#' @slot p.lm.form List of dsign matrices based on formula
#' @slot phi.lm.form List of dsign matrices based on formula
#' @slot lam.lm.form List of dsign matrices based on formula
#' @slot MR_settings A list of mark-recapture settings
#' @export
setClass("tmb_list",
         slots = list(
           parameters = 'list',
           TMB = 'list',
           MR_settings = 'list',
           p.lm.form = 'list',
           phi.lm.form = 'list',
           lam.lm.form = 'list',
           plot = 'list'
         )
)

#'TMBmodel class
#' @slot TMB The TMB obj
#' @export
setClass(
  "tmb_model",
  contains = "Model",
  slots = list(
    TMB = "list"
  )
)
