\name{plot-methods}
\alias{plot}
\alias{PlotMenu}

\title{Plot output from an object of the from the \link{GAS} package.}

\description{

  This method provides plot functionalities for the \link{uGASFit}, \link{mGASFit}, \link{uGASSim},
  \link{mGASSim}, \link{uGASFor}, \link{mGASFor}, \link{uGASRoll} and \link{mGASRoll} objects defined
  in the \link{GAS} package.
}


\usage{
  plot(x, y, \ldots)
%
  PlotMenu(x)
}
%
%\S4method{plot}{uGASFit, mGASFit, uGASSim, mGASSim, uGASFor, mGASFor, uGASRoll, mGASRoll} plot(x, which = NULL, \ldots)
\arguments{

\item{x, y}{
objects of class  \link{uGASFit}, \link{mGASFit}, \link{uGASSim},
  \link{mGASSim}, \link{uGASFor}, \link{mGASFor}, \link{uGASRoll}, \link{mGASRoll}.
}

\item{\ldots}{additional arguments, see Details}

}

\details{

\code{plot} accepts the additional argument \code{numeric} argument \code{which}. By default \code{which = NULL}. If \code{which} is provided, \code{plot()} does not show the interactive menu and plot the corresponding option. The available options for each object class is printed by the function \code{PlotMenu(x)}. By default \code{which = NULL}, that is, \code{plot()} display an interactive menu.
}

\value{

Displays a plot of an object of class  \link{uGASFit}, \link{mGASFit}, \link{uGASSim},
  \link{mGASSim}, \link{uGASFor}, \link{mGASFor}, \link{uGASRoll}, \link{mGASRoll}.

}

\author{Leopoldo Catania}
\examples{
\dontrun{
## Plot filtered estimates of a GAS model estimated on the
## Quarterly logarithmic change in percentage points of the Consumer Price Index data set (cpichg)
library("GAS")

data("cpichg")

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity",
                     GASPar = list(location = TRUE, scale = TRUE,
                                   shape = FALSE))

Fit = UniGASFit(GASSpec, cpichg)

plot(Fit, which = 1)
}
}
