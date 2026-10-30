\name{plot-methods}
\alias{plot}
\alias{plot,uGASFit,missing-method}
\alias{plot,mGASFit,missing-method}
\alias{plot,uGASSim,missing-method}
\alias{plot,mGASSim,missing-method}
\alias{plot,uGASFor,missing-method}
\alias{plot,mGASFor,missing-method}
\alias{plot,uGASRoll,missing-method}
\alias{plot,mGASRoll,missing-method}
\alias{PlotMenu}

\title{Plot output from an object of the \link{GAS} package.}

\description{

  This method provides plot functionalities for the \link{uGASFit}, \link{mGASFit}, \link{uGASSim},
  \link{mGASSim}, \link{uGASFor}, \link{mGASFor}, \link{uGASRoll} and \link{mGASRoll} objects defined
  in the \link{GAS} package.
}


\usage{
  \S4method{plot}{uGASFit,missing}(x, which = NULL)
  \S4method{plot}{mGASFit,missing}(x, which = NULL)
  \S4method{plot}{uGASSim,missing}(x, which = NULL)
  \S4method{plot}{mGASSim,missing}(x, which = NULL)
  \S4method{plot}{uGASFor,missing}(x, which = NULL)
  \S4method{plot}{mGASFor,missing}(x, which = NULL)
  \S4method{plot}{uGASRoll,missing}(x, which = NULL)
  \S4method{plot}{mGASRoll,missing}(x, which = NULL)
%
  PlotMenu(x)
}
%
\arguments{

\item{x}{
objects of class  \link{uGASFit}, \link{mGASFit}, \link{uGASSim},
  \link{mGASSim}, \link{uGASFor}, \link{mGASFor}, \link{uGASRoll}, \link{mGASRoll}.
}

\item{which}{Numeric menu option to plot. If \code{NULL}, display the interactive menu. See Details.}

}

\details{

\code{plot} accepts the additional numeric argument \code{which}. By default \code{which = NULL}. If \code{which} is provided, \code{plot()} does not show the interactive menu and plot the corresponding option. The available options for each object class is printed by the function \code{PlotMenu(x)}. By default \code{which = NULL}, that is, \code{plot()} display an interactive menu.
}

\value{

Displays a plot of an object of class  \link{uGASFit}, \link{mGASFit}, \link{uGASSim},
  \link{mGASSim}, \link{uGASFor}, \link{mGASFor}, \link{uGASRoll}, \link{mGASRoll}.

}

\author{Leopoldo Catania}
\examples{

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
