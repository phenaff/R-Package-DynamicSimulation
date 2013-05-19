\name{UnknownVolLogNormal}
\alias{UnknownVolLogNormal}
\title{Stochastic Volatility path generator}
\usage{
  UnknownVolLogNormal(S0, eps, delta.t, params)
}
\arguments{
  \item{S0}{[numeric] initial value}

  \item{eps}{[numeric] matrix of uniform (0,1) deviates}

  \item{delta.t}{[numeric] time step}

  \item{params}{[list] process definition: \describe{
  \item{mu [numeric]}{drift (annual rate)} \item{sigmaMin
  [numeric]}{min standard deviation (annual rate)}
  \item{sigmaMax [numeric]}{max standard deviation (annual
  rate)} }}
}
\value{
  matrix of paths, one path by column
}
\description{
  Log-normal dynamic with constant drift and stochastic
  volatility. volatility is drawn uniformly between
  \code{sigmaMin} and \code{sigmaMax}
}

