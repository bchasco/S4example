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
#' @slot vars Predictor and response variables
#' @slot lm.form List of dsign matrices based on formula
#' @slot p.lm.form List of dsign matrices based on formula
#' @slot phi.lm.form List of dsign matrices based on formula
#' @slot lam.lm.form List of dsign matrices based on formula
#' @slot formula A formula for the model
#' @slot p.frm A formula for the model
#' @slot phi.frm A formula for the model
#' @slot phi.tmp.frm A formula for the model
#' @slot lam.tmp.frm A formula for the model
#' @slot p.tmp.frm A formula for the model
#' @slot lam.frm A formula for the model
#' @slot phi.list A list of predictor matrices for phi
#' @slot lam.list A list of predictor matrices for phi
#' @slot MR_settings A list of mark-recapture settings
#' @slot plot A list of plotting objects
#' @export
setClass("tmb_list",
         slots = list(
           raw_data = "data.frame",
           parameters = 'list',
           vars = 'list',
           TMB = 'list',
           MR_settings = 'list',
           formula = 'character',
           p.lm.form = 'list',
           phi.lm.form = 'list',
           lam.lm.form = 'list',
           lm.form = 'list',
           phi.list = "list",
           lam.list = 'list',
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
