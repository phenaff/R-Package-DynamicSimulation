#' Delta hedging simulator
#'
#' This function simulates a dynamic hedging strategy of a derivative
#' or of a portfolio of derivatives, all function of the same underlying asset.
#' @title Dynamic delta hedging
#' @param instruments list of instruments to be hedged
#' @param env a \code{\linkS4class{DataProvider}}
#' @param params list of parameters that define the hedging policy:
#' \describe{
#' \item{dtSim [vector of timeDate]}{simulation dates}
#' \item{dtFirst [timeDate]}{first simulation date}
#' \item{dtLast [timeDate]}{last simulation date}
#' \item{nbSteps [numeric]}{number of simulation steps}
#' \item{transaction.cost [numeric]}{proportional transaction cost for the underlying asset}
#' }
#' @param trace output debuging information?
#' @return a list with the following items:
#' \describe{
#' \item{wealth [matrix]}{residual wealth by time step and scenario}
#' \item{portfolio [matrix]}{value of hedge portfolio by time step and scenario}
#' \item{bond [matrix]}{quantity of zero-coupon bond held in the hedge portfolio by time step and scenario}
#' \item{price [matrix]}{price of derivative portfolio, by time step and scenario}
#' \item{stock [matrix]}{quantity of underlying asset in the hedge portfolio, by time step and scenario}
#' \item{description [string]}{description of hedging strategy}
#' }
#'
#' @examples
#' library(fInstrument)
#'
#' dtExpiry <- mytDate('01jan2011')
#' dtStart <- mytDate('01jan2010')
#' nbSteps <- 100;
#' nbPaths <- 500;
#'
#' dtSim <- timeSequence(dtStart, dtExpiry, length.out=nbSteps+1)
#' horizon <- tDiff(dtExpiry, dtStart)
#' delta.t <- horizon/nbSteps
#'
#' sigma <- .3
#'
#' underlying <- 'IBM'
#' K<-100
#' # define derivative
#' a <- fInstrumentFactory("vanilla", quantity=1,
#'                   params=list(cp='c', strike=K,
#'                     dtExpiry=dtExpiry, 
#'                     underlying=underlying,
#'                     discountRef='USD.LIBOR', trace=FALSE))
#'
#' # market data in default environment - basic PV calculation
#'
#' base.env <- DataProvider()
#' setData(base.env, underlying, 'Price', dtStart, 100)
#' setData(base.env, underlying, 'DivYield', dtStart, .02)
#' setData(base.env, underlying, 'ATMVol', dtStart, sigma)
#' setData(base.env, underlying, 'discountRef', dtStart, 'USD.LIBOR')
#' setData(base.env, 'USD.LIBOR', 'Yield', dtStart, .02)
#'
#' getValue(a, 'Price', dtStart, base.env)
#'
#' # price paths 
#' tSpot <- pathSimulator(dtSim = dtSim, nbPaths=nbPaths, 
#'     innovations.gen=sobolInnovations, path.gen=logNormal, 
#'     path.param = list(mu=0, sigma=sigma), S0=100, antithetic = FALSE, 
#'     standardization = TRUE, trace = FALSE)
#'
#' # derived environment for scenario analysis
#' sce.env <- DataProvider(parent=base.env)
#' setData(sce.env, underlying, 'Price',
#'         time(tSpot), as.matrix(tSpot))
#'
#' # simulate a delta-hedge strategy along each path

#' assets = list(a)
#' res <- deltaHedge(assets, sce.env,
#'                   params=list(dtSim=time(tSpot),
#'                   transaction.cost=0), trace=FALSE)
#'
#' x11()
#' hist(tail(res$wealth,1), 50, xlab="wealth",
#'      main=paste("distribution of wealth at expiry ", attr(a,'desc')))
#' @export

deltaHedge <- function(instruments, env, params, trace=FALSE) {

  
dtSim <- params[['dtSim']]
if(is.null(dtSim)) {
  stop('dtSim argument is missing in params list')
  }
transaction.cost <- params[['transaction.cost']]
if(is.null(transaction.cost)) {
  transaction.cost <- 0
}

if(trace) {
  print(params)
}

if(is.null(dtSim)) {
  nbSteps <- params[['nbSteps']]
  dtStart <- params[['dtStart']]
  dtEnd <- params[['dtEnd']]
  dtSim <- seq(from=dtStart, to=dtEnd, length.out=nbSteps)
}
else {
  nbSteps <- length(dtSim)
}

if(trace) {
  print(dtSim)
  print(paste('nbSteps: ', nbSteps))
}
  
# initialize 

tPrice <- 0
tDelta <- 0
for (a in instruments) {
    tPrice <- tPrice + getValue(a, 'Price', dtSim[1], env)
    tDelta <- tDelta + getValue(a, 'Delta', dtSim[1], env)
}

for(i in seq(2, nbSteps)) {
    p.tmp <- 0
    d.tmp <- 0
    for(a in instruments) {
      p.tmp <- p.tmp + getValue(a, 'Price', dtSim[i], env)
      d.tmp <- d.tmp + getValue(a, 'Delta', dtSim[i], env)
    }
    if(trace) {
      print(paste('i: ', i))
      print(p.tmp)
      print(d.tmp)
    }
    tPrice <- rbind(tPrice, p.tmp)
    tDelta <- rbind(tDelta, d.tmp)
}

if(trace) {
  cat("price...\n")
  print(tPrice)
  cat("delta...\n")
  print(tDelta)
}

# value of hedge portfolio
tV <- matrix(nrow=nbSteps, ncol=ncol(tPrice))
# bond position
tAlpha <- matrix(nrow=nbSteps, ncol=ncol(tPrice))
# spot process
tSpot <- matrix(nrow=nbSteps, ncol=ncol(tPrice))

# discount curve name
discountRef <- instruments[[1]]@params$discountRef
Underlying <- instruments[[1]]@params$underlying

# discount rate
iRate <- getData(env, discountRef, 'Yield', dtSim[1])

# form initial portfolio
tV[1,] <- tPrice[1,]
tSpot[1,] <- getData(env, Underlying, 'Price', dtSim[1])
df <- exp(-iRate*tDiff(dtSim[1], dtSim[nbSteps]))
tAlpha[1,] <- (tV[1,] - tDelta[1,]*tSpot[1,]) / df

for(i in seq(2, nbSteps)) {
  iRate <- getData(env, discountRef, 'Yield', dtSim[i])
  df <- exp(-iRate*tDiff(dtSim[i], dtSim[nbSteps]))
  tSpot[i,] <- getData(env, Underlying, 'Price', dtSim[i])
  tV[i,] <- as.numeric(tAlpha[i-1,]) * df + as.numeric(tDelta[i-1,]) * tSpot[i,]
  # reduce portfolio value by transaction cost
  tV[i,] <- tV[i,] - abs(as.numeric(tDelta[i-1,])-as.numeric(tDelta[i,]))*tSpot[i,]*transaction.cost
  tAlpha[i,] <- (tV[i,] - tDelta[i,] * tSpot[i,])/df
}

  if(trace) {
    cat('tV...\n')
    print(tV)
    cat('tPrice...\n')
    print(tPrice)
  }
res <- list("wealth"= tV-tPrice, "portfolio"=tV, "bond"=tAlpha, "price"=tPrice, "stock"=tDelta, "description"='Delta Hedge', 'spot'=tSpot)

res 
}

#' Detailed Report
#'
#' Builds a detailed report by time step as an xtable, with 6 columns:
#' \enumerate{
#' \item time step
#' \item stock price
#' \item delta
#' \item option value
#' \item bond position
#' \item portfolio value
#' }
#'
#' @param iScenario number of scenario to be displayed
#' @param res result from a dynamic hedging simulation function, such as \code{\link{deltaHedge}}
#' @param dec.digits number of digits to be displayed: one number or an array of length 5.
#' @export

makeTable <- function(iScenario, res, dec.digits=2) {
tSpot <- res$spot
nbSteps <- dim(tSpot)[1]
z <- matrix(nrow=nbSteps, ncol=6)
z[,1] <- seq(1,nbSteps)
z[,2] <- tSpot[,iScenario]
# stock in hedge portfolio
z[,3] <- res$stock[,iScenario]
z[,4] <- res$price[,iScenario]
z[,5] <- res$bond[,iScenario]
z[,6] <- res$portfolio[,iScenario]
colnames(z) <- c("time", "stock price", "delta", "option", "bond pos", "hedge port.")
xt <- xtable(z, label=paste('tab:sim-',iScenario,sep=''), caption='Delta hedging simulation')
digits(xt) = dec.digits
xt
}
