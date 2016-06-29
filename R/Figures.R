PlotMenu<-function(object){
  if(is(object,"uGASFit")){
    vplotMenu = c("Filtered Parameters", "Conditional Moments", "Data")
  }
  if(is(object,"mGASFit")){
    vplotMenu = c("Filtered Parameters", "Conditional Moments", "Data")
  }
  if(is(object,"uGASSim")){
    vplotMenu = c("Filtered Parameters", "Conditional Moments", "Simulated Data")
  }
  if(is(object,"mGASSim")){
    vplotMenu = c("Filtered Parameters", "Conditional Moments", "Simulated Data")
  }
  if(is(object,"uGASFor")){
    Roll = object@Info$Roll
    if(!Roll)  vplotMenu = c("Parameters Forecast", "Parameters Forecast + Filtered Values", "Moments Forecast + In Sample Moments")
    if(Roll)   vplotMenu = c("Parameters Forecast", "Forecast vs Realized", "Moments")
  }
  return(vplotMenu)
}

PlotMultipleSeries<-function(mTheta,iK,iT,vDates){
  if(iK<=5){
    layout(matrix(1:iK,iK,1)
           ,heights=c(rep(2,iK-1),2.5))
    for(i in 1:(iK)){
      if(i==1)         par(mar = c(0,4,0.1,2))
      if(i!=1 & i!=iK) par(mar = c(0,4,0  ,2))
      if(i==iK)        par(mar = c(3,4,0,  2))


        vLim = c(min(mTheta[,i]),max(mTheta[,i]))


      plot(vDates,mTheta[,i],type = "n", xaxt="n", xlab="",ylab="", las=1,
           ylim = vLim)
      grid(nx = 10, ny = 10, col = "gray", lty = "dotted")

      lines(vDates, mTheta[,i], col = "black")
      axis(4,at = mean(vLim), labels = colnames(mTheta)[i],tick = F,padj = -1)
    }
    if(!is(vDates,"integer")){
      axis.Date(1, at=seq(min(vDates), max(vDates), "year"))
      axis.Date(1, at=seq(min( vDates), max(vDates), "quarter"),
                labels = FALSE, tcl = -0.2)
    }else{
      foo = c(1,seq(0,iT,ceiling(iT/20))[-1])
      axis(1,at = foo, labels = (1:iT)[foo])
    }
  }else{
    nPlot = ceiling(iK/10)
    plotSeq = seq(1,iK+1,5)
    Start   = 1.0
    PlotType2 = ""

    for(j in 1:nPlot){
      if(PlotType2!="0"){

        layout(matrix(1:10,5,2)
               ,heights=c(rep(2,4),2.5,rep(2,4),2.5))

        for(i in Start:(Start + 9)){
          if(i<=iK){
            if(any(i==plotSeq))                     par(mar = c(0,4,0.1,2))
            if(all(i!=plotSeq) & all(i!=plotSeq-1)) par(mar = c(0,4,0  ,2))
            if(any(i==plotSeq-1))                   par(mar = c(3,4,0,  2))

            vLim = c(min(mTheta[,i]),max(mTheta[,i]))

            plot(vDates,mTheta[,i],type = "n", xaxt="n", xlab="",ylab="", las=1,
                 ylim = vLim)

            grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
            lines(vDates , mTheta[,i], col = "black")
            axis(4,at = mean(vLim), labels = colnames(mTheta)[i],tick = F,padj = -1)

            if(any(i==plotSeq-1) | (i == iK)){

              if(!is(vDates,"integer")){
                axis.Date(1, at=seq(min(vDates), max(vDates), "year"))
                axis.Date(1, at=seq(min( vDates), max(vDates), "quarter"),
                          labels = FALSE, tcl = -0.2)
              }else{
                foo = c(1,seq(0,iT,ceiling(iT/20))[-1])
                axis(1,at = foo, labels = (1:iT)[foo])
              }

            }
          }
        }
        Start = Start+10
        if(j<nPlot) PlotType2=readline("Print enter for next figures or 0 to exit\n:")

      }
    }

  }

}

PlotSingleSeries<-function(vTheta,iT,vDates){

  vLim = c(min(vTheta),max(vTheta))

  layout(matrix(1,1,1)); par(mar = c(3,4,1,2))

  plot(vDates,vTheta,type = "n", xaxt="n", xlab="",ylab="", las=1,
       ylim = vLim)
  grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
  lines(vDates, vTheta, col = "black")
  axis(4,at = mean(vLim), labels = "vY",tick = F,padj = -1,las = 1)

  if(!is(vDates,"integer")){
    axis.Date(1, at=seq(min(vDates), max(vDates), "year"))
    axis.Date(1, at=seq(min( vDates), max(vDates), "quarter"),
              labels = FALSE, tcl = -0.2)
  }else{
    foo = c(1,seq(0,iT,ceiling(iT/20))[-1])
    axis(1,at = foo, labels = (1:iT)[foo])
  }

}

PlotForecastVsRealized_Univ<-function(mRealVsForecast,vDates_Out){

  layout(matrix(1:1,1,1)
         ,heights=c(4))
  par(mar = c(3,4,0.1,2))

  vLim = c(min(mRealVsForecast),max(mRealVsForecast))

  plot(vDates_Out,mRealVsForecast[,1],type = "n", xaxt="n", xlab="",ylab="", las=1,
                  ylim = vLim)
  grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
  lines(vDates, mRealVsForecast[,2], col = "black")
  lines(vDates, mRealVsForecast[,1], col = "red")

  legend("topright",legend = c("Realized","Predicted"), col = c("black","red"), lty = c(1,1))

  if(!is(vDates_Out,"integer")){
    axis.Date(1, at=seq(min(vDates_Out), max(vDates_Out), "year"))
    axis.Date(1, at=seq(min( vDates_Out), max(vDates_Out),  "quarter"),
              labels = FALSE, tcl = -0.2)
  }else{
    foo = c(1,seq(0,iT,ceiling(iT/20))[-1])
    axis(1,at = foo, labels = (1:iT)[foo])
  }
}

# PlotMultipleSeries_Bands<-function(mTheta,iK,iT,vDates, cBands){
#   if(iK<=5){
#     layout(matrix(1:iK,iK,1)
#            ,heights=c(rep(2,iK-1),2.5))
#     for(i in 1:(iK)){
#       if(i==1)         par(mar = c(0,4,0.1,2))
#       if(i!=1 & i!=iK) par(mar = c(0,4,0  ,2))
#       if(i==iK)        par(mar = c(3,4,0,  2))
#
#
#       vLim = c(min(cBands[,,i]),max(cBands[,,i]))
#
#
#       plot(vDates,mTheta[,i],type = "n", xaxt="n", xlab="",ylab="", las=1,
#            ylim = vLim)
#       grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
#
#       lines(vDates, mTheta[,i], col = "black")
#       axis(4,at = mean(vLim), labels = colnames(mTheta)[i],tick = F,padj = -1)
#     }
#     if(!is(vDates,"integer")){
#       axis.Date(1, at=seq(min(vDates), max(vDates), "year"))
#       axis.Date(1, at=seq(min( vDates), max(vDates), "quarter"),
#                 labels = FALSE, tcl = -0.2)
#     }else{
#       foo = c(1,seq(0,iT,ceiling(iT/20))[-1])
#       axis(1,at = foo, labels = (1:iT)[foo])
#     }
#   }else{
#     nPlot = ceiling(iK/10)
#     plotSeq = seq(1,iK+1,5)
#     Start   = 1.0
#     PlotType2 = ""
#
#     for(j in 1:nPlot){
#       if(PlotType2!="0"){
#
#         layout(matrix(1:10,5,2)
#                ,heights=c(rep(2,4),2.5,rep(2,4),2.5))
#
#         for(i in Start:(Start + 9)){
#           if(i<=iK){
#             if(any(i==plotSeq))                     par(mar = c(0,4,0.1,2))
#             if(all(i!=plotSeq) & all(i!=plotSeq-1)) par(mar = c(0,4,0  ,2))
#             if(any(i==plotSeq-1))                   par(mar = c(3,4,0,  2))
#
#             vLim = c(min(mTheta[,i]),max(mTheta[,i]))
#
#             plot(vDates,mTheta[,i],type = "n", xaxt="n", xlab="",ylab="", las=1,
#                  ylim = vLim)
#
#             grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
#             lines(vDates , mTheta[,i], col = "black")
#             axis(4,at = mean(vLim), labels = colnames(mTheta)[i],tick = F,padj = -1)
#
#             if(any(i==plotSeq-1) | (i == iK)){
#
#               if(!is(vDates,"integer")){
#                 axis.Date(1, at=seq(min(vDates), max(vDates), "year"))
#                 axis.Date(1, at=seq(min( vDates), max(vDates), "quarter"),
#                           labels = FALSE, tcl = -0.2)
#               }else{
#                 foo = c(1,seq(0,iT,ceiling(iT/20))[-1])
#                 axis(1,at = foo, labels = (1:iT)[foo])
#               }
#
#             }
#           }
#         }
#         Start = Start+10
#         if(j<nPlot) PlotType2=readline("Print enter for next figures or 0 to exit\n:")
#
#       }
#     }
#
#   }
#
# }
