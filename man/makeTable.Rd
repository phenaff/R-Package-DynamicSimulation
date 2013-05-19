\name{makeTable}
\alias{makeTable}
\title{Dynamic delta hedging}
\usage{
  makeTable(iScenario, res)
}
\arguments{
  \item{iScenario}{number of scenario to be displayed}

  \item{res}{result from a dynamic hedging simulation
  function, such as \code{\link{deltaHedge}}}
}
\value{
  an xtable with 6 columns: \describe{ \item{time}{time
  step} \item{stock price}{the price of the underlying
  asset} \item{delta}{the delta of the option}
  \item{option}{option value} \item{bond pos}{zero-coupon
  bond position} \item{hedge port.}{net value of hedge
  portfolio} }
}
\description{
  Delta hedging simulator
}
\details{
  This function simulates a dynamic hedging strategy of a
  derivative or of a portfolio of derivatives, all function
  of the same underlying asset.
}

