## Test of delta hedging strategy on Vanilla options

library(fInstrument)
library(DynamicSimulation)
library(tmpHenaff)

doplot <- T

dtExpiry <- mytDate('01jan2011')
dtStart <- mytDate('01jan2010')
nbSteps <- 250;
nbPaths <- 1000;

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

# simulate a delta-hedge strategy along each path

assets = list(a)
res <- deltaHedge(assets, sce.env,
                  params=list(dtSim=time(tSpot),
                  transaction.cost=0), trace=F)

if(doplot) {
# distribution of wealth at horizon
x11()
hist(tail(res$wealth,1), 50, xlab="wealth", main=paste("distribution of wealth at expiry ", a@desc))
}

expiry.error <- tail(res$wealth,1)
print(paste('Mean error: ', round(mean(expiry.error),3), ' sd: ', round(sd(as.vector(expiry.error)),3)))
