setClass("uGASFit",representation(ModelInfo="list",GASDyn="list",Estimates="list", Data = "list"))
setClass("mGASFit",representation(ModelInfo="list",GASDyn="list",Estimates="list", Data = "list"))
setClass("uGASSim",representation(ModelInfo="list",GASDyn="list", Data = "list"))
setClass("mGASSim",representation(ModelInfo="list",GASDyn="list", Data = "list"))
setClass("uGASSpec",representation(Spec="list"))
setClass("mGASSpec",representation(Spec="list"))
setClass("uGASFor",representation(Forecast = "list", Bands = "array", Draws = "matrix",
                                  Info = "list", Data = "list"))
setClass("mGASFor",representation(Forecast = "list", Bands = "array", Draws = "matrix",
                                  Info = "list", Data = "list"))
setClass("uGASRoll",representation(Forecast = "list", Info = "list", Data = "list"))

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

setMethod("show", "uGASFor",
          function(object) {

            iH   = object@Info$iH

            Dist          = getDist(object)
            ScalingType   = getScalingType(object)
            ParNames      = FullNamesUni(Dist)
            iK            = NumberParameters(Dist)
            Roll          = object@Info$Roll
            PointForecast = getForecast(object)

            if(Roll){
              PointForecast = cbind(PointForecast, realized = object@Data$vOut)
            }

            cat(paste("\n------------------------------------------"))
            cat(paste("\n-        Univariate GAS Forecast         -"))
            cat(paste("\n------------------------------------------"))
            cat("\n\nModel Specification")
            cat(paste("\nConditional distribution : ",Dist))
            cat(paste("\nScore scaling type : ",ScalingType))
            cat(paste("\nHorizon : ",iH))
            cat(paste("\nRolling forecast : ",Roll))
            #
            cat(paste("\n------------------------------------------"))
            cat(paste("\nParameters forecast:\n"))
            print(PointForecast)
          }
)

setMethod("show", "uGASRoll",
          function(object) {

            ForecastLength = object@Info$ForecastLength

            Dist          = getDist(object)
            ScalingType   = getScalingType(object)
            ParNames      = FullNamesUni(Dist)
            iK            = NumberParameters(Dist)
            PointForecast = getForecast(object)
            elapsedTime   = object@Info$object

            cat(paste("\n------------------------------------------"))
            cat(paste("\n-    Univariate GAS Rolling Forecast     -"))
            cat(paste("\n------------------------------------------"))
            cat("\n\nModel Specification")
            cat(paste("\nConditional distribution : ",Dist))
            cat(paste("\nScore scaling type : ",ScalingType))
            #
            cat(paste("\n------------------------------------------"))
            cat(paste("\nParameters forecast:\n"))
            print(PointForecast)
            cat(paste("\n------------------------------------------"))
            cat(paste("\n\nElapsed time\t:",round(as.double(elapsedTime,units = "mins"),2),"mins"))
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

setMethod("plot", signature(x='uGASFor',y='missing'),
          function(x,...) {
            iK = x@Info$iK
            vY = x@Data$vY
            iH = x@Info$iH
            iT = length(vY)

            Roll = x@Info$Roll
            vOut = x@Data$vOut

            FilteredParameters = x@Data$FilteredParameters
            FilteredParameters = FilteredParameters[-nrow(FilteredParameters),] #remove one step ahead forecast
            ParametersForecast = getForecast(x)
            cBands             = x@Bands
            Moments            = getMoments(x)

            if(is(vY,"xts")) {
              vDates_In = as.Date(index(vY))
              if(Roll) {
                vDates_Out = as.Date(index(vOut))
                ParametersForecast = xts(ParametersForecast, vDates_Out)
              }
            }else{
              vDates = 1:iT
              if(Roll) vDates_Out = (iT+1):(iT+iH)
            }

            if(dev.cur() != 1) dev.off()
            PlotType=1
            while(PlotType>0){
              if(Roll){
                cat(paste("Print 1-3 or 0 to exit"))
                PlotType=menu(PlotMenu(x))
                if(PlotType == 1) PlotMultipleSeries(ParametersForecast,iK,iH,vDates_Out)
                if(PlotType == 2) {
                  Mu              = Moments[,1]
                  mRealVsForecast = cbind(Mu, vOut)
                  PlotForecastVsRealized_Univ(mRealVsForecast,vDates_Out)
                }
                if(PlotType == 3) {
                  PlotMultipleSeries(Moments,4,iH,vDates_Out)
                }
              }
            }
          }
)



getFilteredParameters = function(object)
{
  UseMethod("getFilteredParameters")
}
.getFilteredParameters<-function(object){
  if(is(object,"uGASFit" ) | is(object,"mGASFit" )) mTheta = object@GASDyn$mTheta
  if(is(object,"uGASSim" ) | is(object,"mGASSim" )) mTheta = object@GASDyn$mTheta

  mTheta = t(mTheta)
  parNames = getParNames(object)
  colnames(mTheta) = parNames

  return(mTheta)
}
setMethod("getFilteredParameters", signature(object = "uGASFit"), .getFilteredParameters)
setMethod("getFilteredParameters", signature(object = "mGASFit"), .getFilteredParameters)
setMethod("getFilteredParameters", signature(object = "uGASSim"), .getFilteredParameters)
setMethod("getFilteredParameters", signature(object = "mGASSim"), .getFilteredParameters)

getObs = function(object)
{
  UseMethod("getObs")
}
.getObs<-function(object){
  if(is(object,"uGASFit" )) Data = object@Data$vY
  if(is(object,"mGASFit" )) Data = object@Data$mY

  if(is(object,"uGASSim" )) Data = object@Data$vY
  if(is(object,"mGASSim" )) Data = object@Data$mY

  if(is(object,"uGASFor" )) Data = object@Data$vY

  if(is(object,"uGASRoll" )) Data = object@Data$vY

  return(Data)
}
setMethod("getObs", signature(object = "uGASFit"), .getObs)
setMethod("getObs", signature(object = "mGASFit"), .getObs)
setMethod("getObs", signature(object = "uGASSim"), .getObs)
setMethod("getObs", signature(object = "mGASSim"), .getObs)
setMethod("getObs", signature(object = "uGASFor"), .getObs)
setMethod("getObs", signature(object = "uGASRoll"), .getObs)


getMoments = function(object)
{
  UseMethod("getMoments")
}
.getMoments<-function(object){
  if(is(object,"uGASFit" )) Moments = object@Estimates$Moments
  if(is(object,"mGASFit" )) Moments = NULL

  if(is(object,"uGASSim" )) Moments = object@Data$Moments
  if(is(object,"mGASSim" )) Moments = NULL

  if(is(object,"uGASFor" )) Moments = object@Forecast$Moments
  if(is(object,"mGASFor" )) Moments = NULL

  if(is(object,"uGASRoll" )) Moments = object@Forecast$Moments
  if(is(object,"mGASRoll" )) Moments = NULL

  if(!is.null(Moments)) colnames(Moments) = paste("M",1:4,sep = "")

  return(Moments)
}
setMethod("getMoments", signature(object = "uGASFit"), .getMoments)
setMethod("getMoments", signature(object = "mGASFit"), .getMoments)
setMethod("getMoments", signature(object = "uGASSim"), .getMoments)
setMethod("getMoments", signature(object = "mGASSim"), .getMoments)
setMethod("getMoments", signature(object = "uGASFor"), .getMoments)
setMethod("getMoments", signature(object = "mGASFor"), .getMoments)
setMethod("getMoments", signature(object = "uGASRoll"), .getMoments)
setMethod("getMoments", signature(object = "mGASRoll"), .getMoments)

.getCoef<-function(object){
  if(is(object,"uGASFit" ) | is(object,"mGASFit" )) ans = list(lCoef = object@Estimates$lParList, mCoef = object@Estimates$Inference)
  if(is(object,"uGASSim" ) | is(object,"mGASSim" )) ans = list(vKappa = object@ModelInfo$vKappa,
                                                               mA = object@ModelInfo$mA,
                                                               mB = object@ModelInfo$mB )
  return(ans)

}

setMethod("coef", signature(object = "uGASFit"), .getCoef)
setMethod("coef", signature(object = "mGASFit"), .getCoef)
setMethod("coef", signature(object = "uGASSim"), .getCoef)
setMethod("coef", signature(object = "mGASSim"), .getCoef)

.getQuantile<-function(x, probs = c(0.01, 0.05)){

  if(is(x,"uGASFit" )) mTheta = getFilteredParameters(x)
  if(is(x,"uGASSim" )) mTheta = getFilteredParameters(x)
  if(is(x,"uGASFor" )) mTheta = getForecast(x)
  if(is(x,"uGASRoll" )) mTheta = getForecast(x)


  Dist = getDist(x)

  mQuantile = Quantiles(t(mTheta), Dist, probs)
  colnames(mQuantile) = paste("q.",probs,sep ="")

  return(mQuantile)

}

setMethod("quantile", signature(x = "uGASFit"), .getQuantile)
setMethod("quantile", signature(x = "uGASSim"), .getQuantile)
setMethod("quantile", signature(x = "uGASFor"), .getQuantile)
setMethod("quantile", signature(x = "uGASRoll"), .getQuantile)


pit = function(object)
{
  UseMethod("pit")
}

setMethod("pit", signature(object = "uGASFit"), function(object) object@Estimates$vU)
setMethod("pit", signature(object = "uGASFor"), function(object) object@Forecast$vU)
setMethod("pit", signature(object = "uGASRoll"), function(object) object@Forecast$vU)


LogScore = function(object)
{
  UseMethod("LogScore")
}

setMethod("LogScore", signature(object = "uGASFor"), function(object) object@Forecast$vLS)
setMethod("LogScore", signature(object = "uGASRoll"), function(object) object@Forecast$vLS)


getForecast = function(object)
{
  UseMethod("getForecast")
}
setMethod("getForecast", signature(object = "uGASFor"), function(object) return(object@Forecast$PointForecast))
setMethod("getForecast", signature(object = "mGASFor"), function(object) return(object@Forecast$PointForecast))
setMethod("getForecast", signature(object = "uGASRoll"), function(object) return(object@Forecast$PointForecast))

