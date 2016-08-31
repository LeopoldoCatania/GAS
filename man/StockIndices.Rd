\name{StockIndices}
\docType{data}
\alias{StockIndices}
\title{Data: Daily logarithmic returns in percentage points of the DAX, FTSEMIB and CAC40 from 2007-01-03 to 2016-06-24}
\description{
Daily logarithmic returns in percentage points of the DAX, FTSEMIB and CAC40 from 2007-01-03 to 2016-06-24 obtained at \url{https://finance.yahoo.com/}.
}
\usage{data("StockIndices")}
\format{A \code{matrix} object of dimension 2,445 x 3 containing the daily logarithmic returns in percentage points from 2007-01-03 to 2016-06-24. Missing values are simply removed.}
\references{
\url{https://finance.yahoo.com/}
}
\keyword{datasets}
\examples{
\dontrun{
library("quantmod")

Ticker = c( "^GDAXI", "FTSEMIB.MI", "^FCHI" )

From = "2007-01-01"
To = "2016-06-24"

StockEnv = new.env(has = TRUE)

getSymbols(Ticker, from = From, to = To, env = StockEnv)

mPrices = do.call(cbind, eapply(StockEnv, Ad ))

mRet = diff( log( mPrices ) )

colnames(mRet) = c( "DAX", "FTSEMIB", "CAC40" )

StockIndices = mRet[-1, ]
}
}
