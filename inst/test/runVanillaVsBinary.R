## Comnpare Delta hedge of Vanilla and Binary 

library(fInstrument)
library(DynamicSimulation)
library(tmpHenaff)
library('BB')

doplot <- T

dtExpiry <- myDate('01jan2011')
dtStart <- myDate('01jan2010')
nbSteps <- 100;
nbPaths <- 5000;

dtSim <- seq(as.timeDate(dtStart), as.timeDate(dtExpiry),
             length.out=nbSteps+1)
sigma <- .3

# Vanilla Call
underlying <- 'IBM'
K<-100
a1 <- fInstrumentFactory("vanilla", quantity=1,
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

P <- getValue(a1, 'Price', dtStart, base.env)
D <- getValue(a1, 'Delta', dtStart, base.env)

print(paste('Price, Delta of Vanilla Call:', P, D))

# compute quantity and strike to match price and delta 

foo <- function(x, base.env) {
 q <- x[1]
 K <- x[2]
 a <- fInstrumentFactory("binary", quantity=q,
                  params=list(cp='c', strike=K,
                    dtExpiry=dtExpiry, 
                    underlying=underlying,
                    discountRef='USD.LIBOR', trace=F))

 Pa <- getValue(a, 'Price', dtStart, base.env)
 Da <- getValue(a, 'Delta', dtStart, base.env)

 f <- rep(NA, 2) 
 f[1] <- (P-Pa) 
 f[2] <- (D-Da)
 f
}

ans <- dfsane(par=c(100, 100), fn=foo, base.env=base.env)

q <- ans$par[1]
K <- ans$par[2]

# Binary Call
a2 <- fInstrumentFactory("binary", quantity=q,
                  params=list(cp='c', strike=K,
                    dtExpiry=dtExpiry, 
                    underlying=underlying,
                    discountRef='USD.LIBOR', trace=F))

# plot of price as a function of underlying value 

if(doplot) { 
sce <- t(as.matrix(seq(60, 250, length.out=100)))

# derived data provider with price scenarii
tmp.env <- DataProvider(parent=base.env)
setData(tmp.env, underlying, 'Price', dtStart, sce)

# Compute and plot NPV per scenario for underlying spot price
p1 <- getValue(a1, 'Price', dtStart, tmp.env)
p2 <- getValue(a2, 'Price', dtStart, tmp.env)
plot(sce, p1,main='Vanilla and Binary Prices', type='l', lwd=2, ylab='Price', xlab='Spot', col='green')
lines(sce, p2, type='l', lwd=2, col='red') 
}

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

assets = list(a1)
res1 <- deltaHedge(assets, sce.env,
                  params=list(dtSim=time(tSpot),
                  transaction.cost=0), trace=F)

assets = list(a2)
res2 <- deltaHedge(assets, sce.env,
                  params=list(dtSim=time(tSpot),
                  transaction.cost=0), trace=F)

if(doplot) {
par(mfrow=c(1,2))
# distribution of wealth at horizon
hist(tail(res1$wealth,1), 50, xlab="wealth", main=paste("distribution of wealth at expiry ", a1@desc))
hist(tail(res2$wealth,1), 50, xlab="wealth", main=paste("distribution of wealth at expiry ", a2@desc))
par(mfrow=c(1,1))
}
