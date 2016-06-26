#' FitUGAS Class
#'
#' Class for estimated models belonging to the Univariate Generalized Autoregressive Score family
#'
#' @slot model Contains all the model information
#' @slot data contains all the data in input and output
#' @slot fit contains all the information about fit such as estimates
#' @export
setClass("uGASFit",representation(ModelInfo="list",GASDyn="list",Estimates="list", Data = "list"))

#' FitMGAS Class
#'
#' Class for estimated models belonging to the Multivariate Generalized Autoregressive Score family
#'
#' @slot model Contains all the model information
#' @slot data contains all the data in input and output
#' @slot fit contains all the information about fit such as estimates
#' @export
setClass("mGASFit",representation(ModelInfo="list",GASDyn="list",Estimates="list", Data = "list"))

#' SimUGAS Class
#'
#' Class for simulated Univariate Generalized Autoregressive Score models
#'
#' @slot model Contains all the model information
#' @slot data contains all the data in input and output
#' @slot fit contains all the information about fit such as estimates
#' @export
setClass("uGASSim",representation(ModelInfo="list",GASDyn="list", Data = "list"))

#' SimMGAS Class
#'
#' Class for simulated Multivariate Generalized Autoregressive Score models
#'
#' @slot model Contains all the model information
#' @slot data contains all the data in input and output
#' @slot fit contains all the information about fit such as estimates
#' @export
setClass("mGASSim",representation(ModelInfo="list",GASDyn="list", Data = "list"))

#' specuGAS Class
#'
#' Class for Univariate Generalized Autoregressive Score model specification
#'
#' @slot Spec Contains all the model information
#' @export
setClass("uGASSpec",representation(Spec="list"))

#' SimMGAS Class
#'
#' Class for simulated Multivariate Generalized Autoregressive Score models
#'
#' @slot Spec Contains all the model information
#' @export
setClass("mGASSpec",representation(Spec="list"))

setMethod("show", "uGASFit",
          function(object) {

            Spec = getSpec(object)
            Dist = getDist(object)
            iT   = object@ModelInfo$iT
            iK   = NumberParameters(Dist)
            IC   = getIC(object)

            ParNames    = FullNamesUni(Dist)

            ScalingType = getScalingType(Spec)
            GASPar      = unlist(getGASPar(Spec))
            GASPar      = names(GASPar[GASPar])

            Inference = object@Estimates$Inference

            vTheta_Unc = c(MapParameters_univ(object@Estimates$lParList$vKappa, Dist,iK)); names(vTheta_Unc) = ParNames

            elapsedTime = object@ModelInfo$elapsedTime

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
            cat(paste("\nUnconditional Parameters:\n"))
            print(vTheta_Unc)
            cat(paste("\n------------------------------------------"))
            cat(paste("\nInformation Criteria:\n"))
            print(IC)
            cat(paste("\n------------------------------------------"))
            cat(paste("\n\nElapsed time\t:",round(as.double(elapsedTime,units = "mins"),2),"mins"))
          }
)

setMethod("show", "mGASFit",
          function(object) {

            Spec = getSpec(object)
            iT   = object@ModelInfo$iT
            iN   = object@ModelInfo$iN
            iK   = object@ModelInfo$iK
            IC   = getIC(object)

            Dist        = getDist(Spec)
            ScalingType = getScalingType(Spec)
            GASPar      = unlist(getGASPar(Spec))
            GASPar      = names(GASPar[GASPar])

            ParNames = FullNamesMulti(iN, Dist)

            vTheta_Unc  = c(MapParameters_multi(object@Estimates$lParList$vKappa, Dist,iN,iK)); names(vTheta_Unc) = ParNames
            Inference   = object@Estimates$Inference
            elapsedTime = object@ModelInfo$elapsedTime

            cat(paste("\n------------------------------------------"))
            cat(paste("\n-        Multivariate GAS Fit            -"))
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
            cat(paste("\nUnconditional Parameters:\n"))
            print(vTheta_Unc)
            cat(paste("\n------------------------------------------"))
            cat(paste("\nInformation Criteria:\n"))
            print(IC)
            cat(paste("\n------------------------------------------"))
            cat(paste("\n\nElapsed time\t:",round(as.double(elapsedTime,units = "mins"),2),"mins"))
          }
)

setMethod("show", "uGASSim",
          function(object) {

            iT   = object@ModelInfo$iT

            Dist        = getDist(object)
            ScalingType = getScalingType(object)
            ParNames    = FullNamesUni(Dist)
            iK          = NumberParameters(Dist)

            mA     = object@ModelInfo$mA
            mB     = object@ModelInfo$mA
            vKappa = object@ModelInfo$vKappa ; names(vKappa) = ParNames

            vTheta_Unc = c(MapParameters_univ(vKappa, Dist, iK)) ; names(vTheta_Unc) = ParNames

            cat(paste("\n------------------------------------------"))
            cat(paste("\n-          Univariate GAS Sim            -"))
            cat(paste("\n------------------------------------------"))
            cat("\n\nModel Specification\t:\t")
            cat(paste("\nT = ",iT))
            cat(paste("\nConditional distribution : ",Dist))
            cat(paste("\nScore scaling type : ",ScalingType))
            #
            cat(paste("\n------------------------------------------"))
            cat(paste("\nParameters:\n"))
            cat("vKappa:\n")
            print(vKappa)
            cat("mA:\n")
            print(mA)
            cat("mB:\n")
            print(mB)
            cat(paste("\n------------------------------------------"))
            cat(paste("\nUnconditional Parameters:\n"))
            print(vTheta_Unc)
          }
)

setMethod("show", "mGASSim",
          function(object) {

            iT   = object@ModelInfo$iT
            iN   = object@ModelInfo$iN

            Dist        = getDist(object)
            ScalingType = getScalingType(object)
            ParNames    = FullNamesMulti(iN,Dist)
            iK          = NumberParameters(Dist,iN)

            mA     = object@ModelInfo$mA
            mB     = object@ModelInfo$mA
            vKappa = object@ModelInfo$vKappa; names(vKappa) = ParNames

            vTheta_Unc = c(MapParameters_multi(vKappa, Dist, iN, iK)) ; names(vTheta_Unc) = ParNames

            cat(paste("\n------------------------------------------"))
            cat(paste("\n-          Univariate GAS Sim            -"))
            cat(paste("\n------------------------------------------"))
            cat("\n\nModel Specification\t:\t")
            cat(paste("\nT = ",iT))
            cat(paste("\nN = ",iN))
            cat(paste("\nConditional distribution : ",Dist))
            cat(paste("\nScore scaling type : ",ScalingType))
            #
            cat(paste("\n------------------------------------------"))
            cat(paste("\nParameters:\n"))
            cat("vKappa:\n")
            print(vKappa)
            cat("mA:\n")
            print(mA)
            cat("mB:\n")
            print(mB)
            cat(paste("\n------------------------------------------"))
            cat(paste("\nUnconditional Parameters:\n"))
            print(vTheta_Unc)
          }
)

setMethod("plot", signature(x='uGASFit',y='missing'),
          function(x,...) {
            iK = x@ModelInfo$iK
            iT = x@ModelInfo$iT
            vY = x@Data$vY

            if(is(vY,"xts")) {vDates = index(vY)}else{vDates = 1:length(vY)}
            if(dev.cur() != 1) dev.off()
            PlotType=1
            while(PlotType>0){

              cat(paste("Print 1-3 or 0 to exit"))
              PlotType=menu(PlotMenu(x))
              if(PlotType==1) series2plot = getFilteredParameters(x)[1:iT,]
              if(PlotType==3) series2plot = vY

              if(PlotType==1) PlotMultipleSeries(series2plot,iK,iT,vDates)
              if(PlotType==3) PlotSingleSeries(series2plot,iT,vDates)

            }
          }
)

setMethod("plot", signature(x='mGASFit',y='missing'),
          function(x,...) {
            iK = x@ModelInfo$iK
            iN =x@ModelInfo$iN
            iT = x@ModelInfo$iT
            mY = t(x@Data$mY)

            if(is(mY,"xts")) {vDates = index(mY)}else{vDates = 1:nrow(mY)}
            if(dev.cur() != 1) dev.off()
            PlotType=1
            while(PlotType>0){

              cat(paste("Print 1-3 or 0 to exit"))
              PlotType=menu(PlotMenu(x))
              if(PlotType==1) series2plot = getFilteredParameters(x)[1:iT,]
              if(PlotType==3) series2plot = mY

              if(PlotType==1) PlotMultipleSeries(series2plot,iK,iT,vDates)
              if(PlotType==3) PlotMultipleSeries(series2plot,iN,iT,vDates)

            }
          }
)


setMethod("plot", signature(x='uGASSim',y='missing'),
          function(x,...) {
            iK = x@ModelInfo$iK
            iT = x@ModelInfo$iT
            vY = x@Data$vY

            if(is(vY,"xts")) {vDates = index(vY)}else{vDates = 1:length(vY)}

            if(dev.cur() != 1) dev.off()
            PlotType=1
            while(PlotType>0){

              cat(paste("Print 1-3 or 0 to exit"))
              PlotType=menu(PlotMenu(x))
              if(PlotType==1) series2plot = getFilteredParameters(x)[1:iT,]
              if(PlotType==3) series2plot = vY

              if(PlotType==1) PlotMultipleSeries(series2plot,iK,iT,vDates)
              if(PlotType==3) PlotSingleSeries(series2plot,iT,vDates)


            }
          }
)


setMethod("plot", signature(x='mGASSim',y='missing'),
          function(x,...) {
            iK = x@ModelInfo$iK
            iT = x@ModelInfo$iT
            mY = t(x@Data$mY)

            if(is(mY,"xts")) {vDates = index(mY)}else{vDates = 1:nrow(mY)}

            if(dev.cur() != 1) dev.off()
            PlotType=1
            while(PlotType>0){

              cat(paste("Print 1-3 or 0 to exit"))
              PlotType=menu(PlotMenu(x))
              if(PlotType==1) series2plot = getFilteredParameters(x)[1:iT,]
              if(PlotType==3) series2plot = mY

              if(PlotType==1) PlotMultipleSeries(series2plot,iK,iT,vDates)
              if(PlotType==3) PlotMultipleSeries(series2plot,iN,iT,vDates)

            }
          }
)

