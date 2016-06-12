#' FitUGAS Class
#'
#' Class for estimated models belonging to the Generalized Autoregressive Score family
#'
#' @slot model Contains all the model information
#' @slot data contains all the data in input and output
#' @slot fit contains all the information about fit such as estimates
#' @export
setClass("uGASFit",representation(ModelInfo="list",GASDyn="list",Estimates="list", Data = "list"))

setMethod("show", "uGASFit",
          function(object) {

            Spec = getSpec(object)
            iT   = object@ModelInfo$iT
            elapsedTime = object@ModelInfo$elapsedTime

            Dist        = Spec$Dist
            ScalingType = Spec$ScalingType
            GASPar      = unlist(Spec$GASPar)
            GASPar      = GASPar[GASPar]

            Inference = object@Estimates$Inference

            cat(paste("\n------------------------------------------"))
            cat(paste("\n-          Univariate GAS Fit            -"))
            cat(paste("\n------------------------------------------"))
            cat("\n\nModel Specification\t:\t")
            cat(paste("\nT = ",iT))
            cat(paste("\nConditional distribution : ",Dist))
            cat(paste("\nScore scaling type : ",ScalingType))
            cat(paste("\nTime varying parameres : ", paste(GASPar, collapse = ", ")))
            #
            cat(paste("\n------------------------------------------"))
            cat(paste("\nEstimates:\n"))
            print(Inference)
            cat(paste("\n------------------------------------------"))
            cat(paste("\n\nElapsed time\t:",round(as.double(elapsedTime,units = "mins"),2),"mins"))
          }
)
