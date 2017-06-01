\name{dji30ret}
\docType{data}
\alias{dji30ret}
\title{data: Dow Jones 30 Constituents Closing Value Log Return in percentage points}
\description{
This dataset is taken from the \code{rugarch} package of Ghalanos (2015). Returns are in percentage points.\cr
%

Dow Jones 30 Constituents closing value log returns from 1987-03-16 to
2009-02-03 from Yahoo Finance. Note that AIG was replaced by KFT (Kraft Foods)
on September 22, 2008. This is not reflected in this data set as that would
bring the starting date of the data to 2001.
}
\usage{data("dji30ret")}
\format{A data.frame containing 5,521x30 observations.}
\source{Yahoo Finance}
\references{
Ghalanos A (2015).
"rugarch: Univariate GARCH models."
\url{https://cran.r-project.org/package=rugarch}.}
\keyword{datasets}
