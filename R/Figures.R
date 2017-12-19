PlotMenu <- function(x) {
    if (is(x, "uGASFit")) {
        vplotMenu = c("Filtered Parameters",
                      "Conditional Moments",
                      "Probability Integral Transform",
                      "Data",
                      "Data + Filtered Mean")
    }
    if (is(x, "mGASFit")) {
        vplotMenu = c("Filtered Parameters",
                      "Conditional Moments",
                      "Data")
    }
    if (is(x, "uGASSim")) {
        vplotMenu = c("Filtered Parameters",
                      "Conditional Moments",
                      "Simulated Data")
    }
    if (is(x, "mGASSim")) {
        vplotMenu = c("Filtered Parameters",
                      "Conditional Moments",
                      "Simulated Data")
    }
    if (is(x, "uGASFor")) {
        Roll = x@Info$Roll
        if (!Roll)
            vplotMenu = c("Parameters Forecast",
                          "Parameters Forecast with confidence bands",
                          "Parameters Forecast + Filtered Values",
                          "Parameters Forecast with confidence bands + Filtered Values",
                          "Conditional Moments Forecast",
                          "Conditional Moments Forecast + In Sample Moments")
        if (Roll)
            vplotMenu = c("Parameters Forecast",
                          "Forecast vs Realized",
                          "Conditional Moments",
                          "Log scores")
    }
    if (is(x, "mGASFor")) {
        Roll = x@Info$Roll
        if (!Roll)
            vplotMenu = c("Parameters Forecast",
                          "Parameters Forecast with confidence bands",
                          "Parameters Forecast + Filtered Values",
                          "Parameters Forecast with confidence bands + Filtered Values")
        if (Roll)
            vplotMenu = c("Parameters Forecast",
                          "Forecast vs Realized",
                          "Conditional Moments",
                          "Log scores")
    }
    if (is(x, "uGASRoll")) {
        vplotMenu = c("Parameters Forecast",
                      "Forecast vs Realized",
                      "Conditional Moments",
                      "Probability Integral Transform")
    }
    if (is(x, "mGASRoll")) {
        vplotMenu = c("Parameters Forecast",
                      "Forecast vs Realized",
                      "Conditional Moments")
    }

    return(vplotMenu)
}

PlotMultipleSeries <- function(mTheta, iK, iT, vDates) {
    if (iK <= 5L) {
        layout(matrix(1:iK, iK, 1L), heights = c(rep(2, iK - 1), 2.5))
        for (i in 1:(iK)) {
            if (i == 1)
                par(mar = c(0, 4, 0.1, 2))
            if (i != 1 & i != iK)
                par(mar = c(0, 4, 0, 2))
            if (i == iK)
                par(mar = c(3, 4, 0, 2))

            vLim = range(na.omit(mTheta[, i]))

            plot(vDates, mTheta[, i], type = "n", xaxt = "n", xlab = "", ylab = "", las = 1, ylim = vLim)
            grid(nx = 10, ny = 10, col = "gray", lty = "dotted")

            lines(vDates, mTheta[, i])
            axis(4, at = mean(vLim), labels = colnames(mTheta)[i], tick = F, padj = -1)
        }
        if (!is(vDates, "integer")) {
            DiffTime = vDates[2] - vDates[1]
            axis.Date(1, at = seq(min(vDates), max(vDates), DiffTime * iT/20))
            axis.Date(1, at = seq(min(vDates), max(vDates), DiffTime * iT/5), labels = FALSE, tcl = -0.2)
        } else {
            foo = vDates[c(1L, seq(0L, iT, ceiling((iT)/20))[-1])]
            axis(1, at = foo, labels = foo)
        }
    } else {
        nPlot = ceiling(iK/10)
        plotSeq = seq(1, iK + 1, 5)
        Start = 1L
        PlotType2 = ""

        for (j in 1:nPlot) {
            if (PlotType2 != "0") {

                layout(matrix(1:10, 5, 2), heights = c(rep(2, 4), 2.5, rep(2, 4), 2.5))

                for (i in Start:(Start + 9)) {
                  if (i <= iK) {
                    if (any(i == plotSeq))
                      par(mar = c(0, 4, 0.1, 2))
                    if (all(i != plotSeq) & all(i != plotSeq - 1))
                      par(mar = c(0, 4, 0, 2))
                    if (any(i == plotSeq - 1))
                      par(mar = c(3, 4, 0, 2))

                    vLim = c(min(mTheta[, i]), max(mTheta[, i]))

                    plot(vDates, mTheta[, i], type = "n", xaxt = "n", xlab = "", ylab = "", las = 1,
                      ylim = vLim)

                    grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
                    lines(vDates, mTheta[, i])
                    axis(4L, at = mean(vLim), labels = colnames(mTheta)[i], tick = F, padj = -1)

                    if (any(i == plotSeq - 1L) | (i == iK)) {

                      if (!is(vDates, "integer")) {
                        axis.Date(1L, at = seq(min(vDates), max(vDates), "year"))
                        axis.Date(1L, at = seq(min(vDates), max(vDates), "quarter"), labels = FALSE,
                          tcl = -0.2)
                      } else {
                        foo = vDates[c(1, seq(0, iT, ceiling((iT)/20))[-1])]
                        axis(1L, at = foo, labels = foo)
                      }

                    }
                  }
                }
                Start = Start + 10L
                if (j < nPlot)
                  PlotType2 = readline("Print enter for next figures or 0 to exit\n:")

            }
        }

    }

}

PlotMultipleSeries_wis <- function(mTheta_is, mTheta_os, iK, iH, vDates_os, vDates_is) {

    iS = nrow(mTheta_is)

    if (is.numeric(vDates_is))
        vDates_os = (length(vDates_is) + 1L):(length(vDates_is) + iH)

    vDateFull = c(vDates_is, vDates_os)

    if (iK <= 5L) {
        layout(matrix(1:iK, iK, 1L), heights = c(rep(2, iK - 1), 2.5))
        for (i in 1:(iK)) {
            if (i == 1L)
                par(mar = c(0, 4, 0.1, 2))
            if (i != 1L & i != iK)
                par(mar = c(0, 4, 0, 2))
            if (i == iK)
                par(mar = c(3, 4, 0, 2))

            vLim = c(min(mTheta_os[, i], mTheta_is[, i]), max(mTheta_os[, i], mTheta_is[, i]))

            plot(vDateFull, rep(0, length(vDateFull)), type = "n", xaxt = "n", xlab = "", ylab = "",
                las = 1L, ylim = vLim)
            grid(nx = 10, ny = 10, col = "gray", lty = "dotted")

            lines(vDates_is, mTheta_is[, i])
            lines(vDates_os, mTheta_os[, i], col = "red")

            abline(v = tail(vDates_is, 1L), lty = 2L)
            axis(4L, at = mean(vLim), labels = colnames(mTheta_is)[i], tick = FALSE, padj = -1)
        }
        if (!is(vDateFull, "integer")) {
            axis.Date(1L, at = seq(min(vDateFull), max(vDateFull), by = vDateFull[2L] - vDateFull[1L]))
            axis.Date(1L, at = seq(min(vDateFull), max(vDateFull), by = 0.25 * (vDateFull[2L] - vDateFull[1L])),
                labels = FALSE, tcl = -0.2)
        } else {
            foo = vDateFull[c(1, seq(0, iS + iH, ceiling((iS + iH)/20))[-1])]
            axis(1, at = foo, labels = foo)
        }
    } else {
        nPlot = ceiling(iK/10)
        plotSeq = seq(1L, iK + 1L, 5)
        Start = 1L
        PlotType2 = ""

        for (j in 1:nPlot) {
            if (PlotType2 != "0") {

                layout(matrix(1:10, 5L, 2L), heights = c(rep(2, 4), 2.5, rep(2, 4), 2.5))

                for (i in Start:(Start + 9L)) {
                  if (i <= iK) {
                    if (any(i == plotSeq))
                      par(mar = c(0, 4, 0.1, 2))
                    if (all(i != plotSeq) & all(i != plotSeq - 1L))
                      par(mar = c(0, 4, 0, 2))
                    if (any(i == plotSeq - 1))
                      par(mar = c(3, 4, 0, 2))

                    vLim = c(min(mTheta_os[, i], mTheta_is[, i]), max(mTheta_os[, i], mTheta_is[, i]))

                    plot(vDateFull, rep(0, length(vDateFull)), type = "n", xaxt = "n", xlab = "", ylab = "",
                      las = 1L, ylim = vLim)
                    grid(nx = 10, ny = 10, col = "gray", lty = "dotted")

                    lines(vDates_is, mTheta_is[, i])
                    lines(vDates_os, mTheta_os[, i], col = "red")

                    abline(v = tail(vDates_is, 1L), lty = 2L)
                    axis(4L, at = mean(vLim), labels = colnames(mTheta_is)[i], tick = FALSE, padj = -1)

                    if (any(i == plotSeq - 1L) | (i == iK)) {

                      if (!is(vDateFull, "integer")) {
                        axis.Date(1L, at = seq(min(vDateFull), max(vDateFull), by = vDateFull[2L] - vDateFull[1L]))
                        axis.Date(1L, at = seq(min(vDateFull), max(vDateFull), by = 0.25 * (vDateFull[2L] -
                          vDateFull[1L])), labels = FALSE, tcl = -0.2)
                      } else {
                        foo = vDateFull[c(1, seq(0, (iS + iH), ceiling((iS + iH)/20))[-1])]
                        axis(1L, at = foo, labels = foo)
                      }

                    }
                  }
                }
                Start = Start + 10L
                if (j < nPlot)
                  PlotType2 = readline("Print enter for next figures or 0 to exit\n:")

            }
        }

    }

}

PlotSingleSeries <- function(vTheta, iT, vDates) {

    vLim = c(min(vTheta), max(vTheta))

    layout(matrix(1, 1L, 1L))
    par(mar = c(3, 4, 1, 2))

    plot(vDates, vTheta, type = "n", xaxt = "n", xlab = "", ylab = "", las = 1L, ylim = vLim)
    grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
    lines(vDates, vTheta)
    axis(4L, at = mean(vLim), labels = "vY", tick = FALSE, padj = -1, las = 1L)

    if (!is(vDates, "integer")) {
        axis.Date(1L, at = seq(min(vDates), max(vDates), "year"))
        axis.Date(1L, at = seq(min(vDates), max(vDates), "quarter"), labels = FALSE, tcl = -0.2)
    } else {
        foo = vDates[c(1, seq(0, iT, ceiling((iT)/20))[-1])]
        axis(1L, at = foo, labels = foo)
    }

}

PlotForecastVsRealized_Univ <- function(mRealVsForecast, vDates_os, object) {

    layout(matrix(1:1, 1, 1), heights = c(4))
    par(mar = c(3, 4, 0.1, 2))

    vLim = c(min(mRealVsForecast), max(mRealVsForecast))

    iH = length(vDates_os)

    plot(vDates_os, mRealVsForecast[, 1L], type = "n", xaxt = "n", xlab = "", ylab = "", las = 1, ylim = vLim)
    grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
    lines(vDates_os, mRealVsForecast[, 2L])
    lines(vDates_os, mRealVsForecast[, 1L], col = "red")

    if (is(object, "uGASFor") | is(object, "uGASRoll")) {
        legend("topright", legend = c("Realized", "Predicted"), col = c("black", "red"), lty = c(1L,
            1L))
    }
    if (is(object, "uGASFit")) {
        legend("topright", legend = c("Realized", "Filtered"), col = c("black", "red"), lty = c(1L, 1L))
    }

    if (!is(vDates_os, "integer")) {
        axis.Date(1L, at = seq(min(vDates_os), max(vDates_os), vDates_os[2L] - vDates_os[1L]))
        axis.Date(1L, at = seq(min(vDates_os), max(vDates_os), 0.25 * (vDates_os[2L] - vDates_os[1L])),
            labels = FALSE, tcl = -0.2)
    } else {
        foo = vDates_os[c(1, seq(0, iH, ceiling((iH)/20))[-1])]
        axis(1, at = foo, labels = foo)
    }
}

PlotForecastVsRealized_Multi <- function(mReal, mForcasted, iN, vDates_os, object) {

    iH = nrow(mForcasted)

    if (iN <= 5L) {
        layout(matrix(1:iN, iN, 1L), heights = c(rep(2, iN - 1), 2.5))

        for (i in 1:(iN)) {
            if (i == 1L)
                par(mar = c(0, 4, 0.1, 2))
            if (i != 1L & i != iN)
                par(mar = c(0, 4, 0, 2))
            if (i == iN)
                par(mar = c(3, 4, 0, 2))

            vLim = c(min(mReal[, i], mForcasted[, i]), max(mReal[, i], mForcasted[, i]))

            plot(vDates_os, mReal[, i], type = "n", xaxt = "n", xlab = "", ylab = "", las = 1L, ylim = vLim)

            grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
            lines(vDates_os, mReal[, i])
            lines(vDates_os, mForcasted[, i], col = "red")

            axis(4L, at = mean(vLim), labels = colnames(mReal)[i], tick = F, padj = -1)

            if (i == 1L) {
                if (is(object, "mGASFor") | is(object, "mGASRoll"))
                  legend("topright", legend = c("Realized", "Predicted"), col = c("black", "red"), lty = c(1L,
                    1L))
                if (is(object, "mGASFit"))
                  legend("topright", legend = c("Realized", "Filtered"), col = c("black", "red"), lty = c(1L,
                    1L))
            }
        }
        if (!is(vDates_os, "integer")) {
            axis.Date(1L, at = seq(min(vDates_os), max(vDates_os), "year"))
            axis.Date(1L, at = seq(min(vDates_os), max(vDates_os), "quarter"), labels = FALSE, tcl = -0.2)
        } else {
            foo = vDates_os[c(1, seq(0, iH, ceiling((iH)/20))[-1])]
            axis(1, at = foo, labels = foo)
        }
    } else {

        nPlot = ceiling(iN/10)
        plotSeq = seq(1, iN + 1, 5)
        Start = 1L
        PlotType2 = ""

        for (j in 1:nPlot) {
            if (PlotType2 != "0") {

                layout(matrix(1:10, 5L, 2L), heights = c(rep(2, 4), 2.5, rep(2, 4), 2.5))
                for (i in Start:(Start + 9L)) {
                  if (i <= iN) {
                    if (any(i == plotSeq))
                      par(mar = c(0, 4, 0.1, 2))
                    if (all(i != plotSeq) & all(i != plotSeq - 1L))
                      par(mar = c(0, 4, 0, 2))
                    if (any(i == plotSeq - 1))
                      par(mar = c(3, 4, 0, 2))

                    vLim = c(min(mReal[, i], mForcasted[, i]), max(mReal[, i], mForcasted[, i]))

                    plot(vDates_os, mReal[, i], type = "n", xaxt = "n", xlab = "", ylab = "", las = 1L,
                      ylim = vLim)

                    grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
                    lines(vDates_os, mReal[, i], col = "black")
                    lines(vDates_os, mForcasted[, i], col = "red")

                    axis(4L, at = mean(vLim), labels = colnames(mReal)[i], tick = FALSE, padj = -1)

                    if (any(i == plotSeq)) {
                      if (is(object, "mGASFor") | is(object, "mGASRoll"))
                        legend("topright", legend = c("Realized", "Predicted"), col = c("black", "red"),
                          lty = c(1L, 1L))
                      if (is(object, "mGASFit"))
                        legend("topright", legend = c("Realized", "Filtered"), col = c("black", "red"),
                          lty = c(1L, 1L))
                    }
                  }
                }
                if (!is(vDates_os, "integer")) {
                  axis.Date(1L, at = seq(min(vDates_os), max(vDates_os), "year"))
                  axis.Date(1L, at = seq(min(vDates_os), max(vDates_os), "quarter"), labels = FALSE,
                    tcl = -0.2)
                } else {
                  foo = vDates_os[c(1L, seq(0L, iH, ceiling((iH)/20))[-1])]
                  axis(1L, at = foo, labels = foo)
                }
                Start = Start + 10L
                if (j < nPlot)
                  PlotType2 = readline("Print enter for next figures or 0 to exit\n:")

            }
        }
    }
}

PlotMultipleSeries_Bands <- function(mTheta, iK, iT, vDates, cBands) {
    iQ = dim(cBands)[2L]
    if (iK <= 5L) {
        layout(matrix(1:iK, iK, 1L), heights = c(rep(2, iK - 1), 2.5))
        for (i in 1:(iK)) {
            if (i == 1L)
                par(mar = c(0, 4, 0.1, 2))
            if (i != 1L & i != iK)
                par(mar = c(0, 4, 0, 2))
            if (i == iK)
                par(mar = c(3, 4, 0, 2))


            vLim = c(min(cBands[, , i]), max(cBands[, , i]))

            plot(vDates, mTheta[, i], type = "n", xaxt = "n", xlab = "", ylab = "", las = 1L, ylim = vLim)
            grid(nx = 10, ny = 10, col = "gray", lty = "dotted")

            lines(vDates, mTheta[, i], col = "red")
            for (q in 1:iQ) {
              lines(vDates, cBands[, q, i], col = "blue")
            }

            axis(4L, at = mean(vLim), labels = colnames(mTheta)[i], tick = FALSE, padj = -1)
        }
        if (!is(vDates, "integer")) {
            axis.Date(1L, at = seq(min(vDates), max(vDates), "year"))
            axis.Date(1L, at = seq(min(vDates), max(vDates), "quarter"), labels = FALSE, tcl = -0.2)
        } else {
            foo = vDates[c(1L, seq(0L, iT, ceiling((iT)/20))[-1])]
            axis(1L, at = foo, labels = foo)
        }
    } else {
        nPlot = ceiling(iK/10)
        plotSeq = seq(1L, iK + 1, 5)
        Start = 1L
        PlotType2 = ""

        for (j in 1:nPlot) {
            if (PlotType2 != "0") {

                layout(matrix(1:10, 5L, 2L), heights = c(rep(2, 4), 2.5, rep(2, 4), 2.5))

                for (i in Start:(Start + 9L)) {
                  if (i <= iK) {
                    if (any(i == plotSeq))
                      par(mar = c(0, 4, 0.1, 2))
                    if (all(i != plotSeq) & all(i != plotSeq - 1))
                      par(mar = c(0, 4, 0, 2))
                    if (any(i == plotSeq - 1))
                      par(mar = c(3, 4, 0, 2))

                    vLim = c(min(cBands[, , i]), max(cBands[, , i]))

                    plot(vDates, mTheta[, i], type = "n", xaxt = "n", xlab = "", ylab = "", las = 1L,
                      ylim = vLim)

                    grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
                    lines(vDates, mTheta[, i])
                    for (q in 1:iQ) {
                      lines(vDates, cBands[, q, i], col = "red")
                    }
                    axis(4L, at = mean(vLim), labels = colnames(mTheta)[i], tick = F, padj = -1)

                    if (any(i == plotSeq - 1L) | (i == iK)) {

                      if (!is(vDates, "integer")) {
                        axis.Date(1L, at = seq(min(vDates), max(vDates), "year"))
                        axis.Date(1L, at = seq(min(vDates), max(vDates), "quarter"), labels = FALSE,
                          tcl = -0.2)
                      } else {
                        foo = vDates[c(1, seq(0, iT, ceiling((iT)/20))[-1])]
                        axis(1L, at = foo, labels = foo)
                      }

                    }
                  }
                }
                Start = Start + 10L
                if (j < nPlot)
                  PlotType2 = readline("Print enter for next figures or 0 to exit\n:")

            }
        }

    }

}

PlotMultipleSeries_Bands_wis <- function(mTheta_is, mTheta_os, iK, iH, vDates_os, vDates_is, cBands) {
    iQ = dim(cBands)[2L]

    iS = nrow(mTheta_is)

    if (is.numeric(vDates_is)) {
        vDates_os = (length(vDates_is) + 1L):(length(vDates_is) + iH)
    }

    vDateFull = c(vDates_is, vDates_os)

    if (iK <= 5L) {
        layout(matrix(1:iK, iK, 1L), heights = c(rep(2, iK - 1), 2.5))
        for (i in 1L:(iK)) {
            if (i == 1L)
                par(mar = c(0, 4, 0.1, 2))
            if (i != 1L & i != iK)
                par(mar = c(0, 4, 0, 2))
            if (i == iK)
                par(mar = c(3, 4, 0, 2))

            vLim = c(min(cBands[, , i], mTheta_is[, i]), max(cBands[, , i], mTheta_is[, i]))

            plot(vDateFull, rep(0L, length(vDateFull)), type = "n", xaxt = "n", xlab = "", ylab = "",
                las = 1L, ylim = vLim)
            grid(nx = 10, ny = 10, col = "gray", lty = "dotted")

            lines(vDates_is, mTheta_is[, i])
            lines(vDates_os, mTheta_os[, i], col = "red")

            for (q in 1:iQ) {
              lines(vDates_os, cBands[, q, i], col = "blue")
            }
            abline(v = tail(vDates_is, 1L), lty = 2L)
            axis(4L, at = mean(vLim), labels = colnames(mTheta_is)[i], tick = FALSE, padj = -1)
        }
        if (!is(vDateFull, "integer")) {
            axis.Date(1L, at = seq(min(vDateFull), max(vDateFull), by = vDateFull[2L] - vDateFull[1L]))
            axis.Date(1L, at = seq(min(vDateFull), max(vDateFull), by = 0.25 * (vDateFull[2L] - vDateFull[1L])),
                labels = FALSE, tcl = -0.2)
        } else {
            foo = vDateFull[c(1L, seq(0L, iH + iS, ceiling((iH + iS)/20))[-1L])]
            axis(1L, at = foo, labels = foo)
        }
    } else {
        nPlot = ceiling(iK/10)
        plotSeq = seq(1L, iK + 1L, 5L)
        Start = 1L
        PlotType2 = ""

        for (j in 1:nPlot) {
            if (PlotType2 != "0") {

                layout(matrix(1:10, 5L, 2L), heights = c(rep(2, 4), 2.5, rep(2, 4), 2.5))

                for (i in Start:(Start + 9L)) {
                  if (i <= iK) {
                    if (any(i == plotSeq)) {
                      par(mar = c(0, 4, 0.1, 2))
                    }
                    if (all(i != plotSeq) & all(i != plotSeq - 1L)) {
                      par(mar = c(0, 4, 0, 2))
                    }
                    if (any(i == plotSeq - 1L)) {
                      par(mar = c(3, 4, 0, 2))
                    }

                    vLim = c(min(cBands[, , i], mTheta_is[, i]), max(cBands[, , i], mTheta_is[, i]))

                    plot(vDateFull, rep(0, length(vDateFull)), type = "n", xaxt = "n", xlab = "", ylab = "",
                      las = 1L, ylim = vLim)
                    grid(nx = 10, ny = 10, col = "gray", lty = "dotted")

                    lines(vDates_is, mTheta_is[, i])
                    lines(vDates_os, mTheta_os[, i], col = "red")

                    for (q in 1:iQ) {
                      lines(vDates_os, cBands[, q, i], col = "blue")
                    }

                    abline(v = tail(vDates_is, 1L), lty = 2L)
                    axis(4L, at = mean(vLim), labels = colnames(mTheta_is)[i], tick = FALSE, padj = -1L)

                    if (any(i == plotSeq - 1L) | (i == iK)) {

                      if (!is(vDateFull, "integer")) {
                        axis.Date(1L, at = seq(min(vDateFull), max(vDateFull), by = vDateFull[2L] - vDateFull[1L]))
                        axis.Date(1L, at = seq(min(vDateFull), max(vDateFull), by = 0.25 * (vDateFull[2L] -
                          vDateFull[1L])), labels = FALSE, tcl = -0.2)
                      } else {
                        foo = vDateFull[c(1L, seq(0L, (iH + iS), ceiling((iH + iS)/20))[-1L])]
                        axis(1L, at = foo, labels = foo)
                      }

                    }
                  }
                }
                Start = Start + 10L
                if (j < nPlot)
                  PlotType2 = readline("Print enter for next figures or 0 to exit\n:")

            }
        }

    }

}

PlotPit <- function(vU, BIN) {

    layout(matrix(1, 1L, 1L), heights = c(4))
    par(mar = c(3, 4, 0.1, 2))

    h = BIN$hist
    n_i = h$counts
    confidence = BIN$confidence

    plot(h, col = "blue", ylim = c(0, max(confidence * 1.2, n_i * 1.2)))
    abline(h = confidence, col = "red", lwd = 2L, xlim = c(0, 1))
}

PlotCovariances <- function(cCov, iN, iT, vDates, vNames) {

    iL = iN * (iN + 1L)/2L

    vNamesFull = character(iL)
    mCov = matrix(0, iT, iL)

    iC = 1L
    for (i in 1:iN) {
        for (j in 1:i) {
            vNamesFull[iC] = paste(vNames[i], vNames[j], sep = "-")
            mCov[, iC] = cCov[i, j, ]
            iC = iC + 1L
        }
    }

    colnames(mCov) = vNamesFull

    PlotMultipleSeries(mCov, iL, iT, vDates)

}
