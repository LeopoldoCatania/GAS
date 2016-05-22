#' FitUGAS Class
#'
#' Class for estimated models belonging to the Generalized Autoregressive Score family
#'
#' @slot model Contains all the model information
#' @slot data contains all the data in input and output
#' @slot fit contains all the information about fit such as estimates
#' @export
setClass("FitUGAS",representation(model="list",data="list",fit="list"))
# #' SpecUGAS Class
# #'
# #' Class for models belonging to the Generalized Autoregressive Score Copula family
# #'
# #' @slot model Contains all the model information
# #' @slot data contains all the data in input and output
# #' @slot fit contains all the information about fit such as estimates
# #' @slot fit contains all the information about fit such as estimates
# #' @import methods
# #' @export
# #' @exportClass
# setClass("SpecUGAS",representation(model="list",data="list",fit="list"))
