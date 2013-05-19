\name{pathSimulator}
\alias{pathSimulator}
\title{Path simulator}
\usage{
  pathSimulator(dtSim = NULL, horizon = NULL,
    delta.t = NULL, nbPaths = NULL,
    innovations.gen = sobolInnovations,
    path.gen = logNormal, path.param, S0 = 100,
    antithetic = TRUE, standardization = FALSE,
    trace = FALSE)
}
\arguments{
  \item{dtSim}{[timeDate] vector of sampling dates}

  \item{horizon}{[numeric] simulation horizon, in fraction
  of years}

  \item{delta.t}{[numeric] time step}

  \item{nbPaths}{[numeric] number of paths to be generated}

  \item{innovations.gen}{[function] generator of uniform
  (0,1) random numbers}

  \item{path.gen}{[function] price dynamic}

  \item{path.param}{[numeric] list of parameters specific
  to \code{path.gen}}

  \item{S0}{[numeric] initial value}

  \item{antithetic}{[boolean] draw antithetic variates?}

  \item{standardization}{[boolean] recenter and normalize
  the output of \code{innovation.gen}?}

  \item{trace}{[boolean] output debugging information?}
}
\value{
  time series of paths, one path per column
}
\description{
  Function for generating saple paths, given a generator of
  uniform (0,1) sequences and the definition of the dynamic
  process
}

