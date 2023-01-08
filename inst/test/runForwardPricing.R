## Test of forward pricing under a number of scenarios

library(fInstrument)
library(DynamicSimulation)
library(tmpHenaff)

doplot <- T

dtExpiry <- myDate('01jan2011')
dtStart <- myDate('01jan2010')
nbSteps <- 10;
nbPaths <- 5;

dtSim <- seq(as.timeDate(dtStart), as.timeDate(dtExpiry),
             length.out=nbSteps+1)
sigma <- .3

underlying <- 'IBM'
K<-100

# define derivative
a <- fInstrumentFactory("vanilla", quantity=1,
                  params=list(cp='c', strike=K,
                    dtExpiry=dtExpiry, 
                    underlying=underlying,
                    discountRef='USD.LIBOR', trace=F))

# market data in default environment - basic PV calculation

base.env <- DataProvider()
setData(base.env, underlying, 'Price', dtStart, 100)
setData(base.env, underlying, 'DivYield', dtStart, .02)
setData(base.env, underlying, 'ATMVol', dtStart, sigma)
setData(base.env, underlying, 'discountRef', dtStart, 'USD.LIBOR')
setData(base.env, 'USD.LIBOR', 'Yield', dtStart, .02)

# value as of calculation date
getValue(a, 'Price', dtStart, base.env)

# price paths 
tSpot <- pathSimulator(dtSim = dtSim, nbPaths=nbPaths, 
    innovations.gen=sobolInnovations, path.gen=logNormal, 
    path.param = list(mu=0, sigma=sigma), S0=100, antithetic = F, 
    standardization = TRUE, trace = F)

# derived environment for scenario analysis
sce.env <- DataProvider(parent=base.env)
setData(sce.env, underlying, 'Price',
        time(tSpot), as.matrix(tSpot))

# vol fetched from base scenario for all dates
# price fetched from child scenario

dtSim = time(tSpot)
for(i in seq_along(dtSim)) {
  p = getData(sce.env, underlying, 'Price', dtSim[i])
  print(p)
  vol = getData(sce.env, underlying, 'ATMVol', dtSim[i])
  print(vol)
}

