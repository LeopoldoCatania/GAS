\name{tqdata}
\docType{data}
\alias{tqdata}
\title{Data from Bien et al (2011).}
\description{
From the readme.bnp.txt file in the JAE Data Archive available at \url{http://qed.econ.queensu.ca/jae/datasets/bien001/}:

The high-frequency data used in the paper come from the Trades and Quotation
(TAQ) database. The data contains time-stamped quotations of Citicorp stock
traded at the NYSE over the period from 20th February to 23rd February 2001.

In the study, 30-second bid and ask quote changes are constructed from the
irregularly-spaced quote data. The study covers observations recorded from
9:35 EST until 16:00 EST.

The data contains 3080 rows and eight columns - in order:

1.	year
2.	month
3.	day
4.	time in number of seconds after the 9:35 EST
5.	best ask quote
6.	best bid quote
7.	30-second change of the ask quote in number of ticks
8.	30-second change of the bid quote in number of ticks.
}
\usage{data("tqdata")}
\format{A \link{data.frame} object containing 3,080 observat. ions.}
\references{
Bien K, Nolte, I, Pohlmeier W (2011). 
"An Inflated Multivariate Integer Count Hurdle Model: An Application to Bid and Ask Quote Dynamics". 
Journal of Applied Econometrics, 26(4), 669-707. 
\doi{10.1002/jae.1122}
}
\keyword{datasets}
