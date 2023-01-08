## Delta Hedging error as a function of hedging frequency
## and transaction cost

library(gregmisc)
library(fInstrument)
library(DynamicSimulation)
library(tmpHenaff)

doplot <- T

dtExpiry <- mytDate('01jan2011')
dtStart <- mytDate('01jan2010')
nbPaths <- 1000;

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

pr <- getValue(a, 'Price', dtStart, base.env)
vega <- getValue(a, 'Vega', dtStart, base.env)
 
print(paste('Initial value of derivative: ' , round(pr,2), 'Vega: ', round(vega,2)))

nb.range <- seq(10, 500, by=50)
mean.error <- 0.0 * nb.range
sd.error <- 0.0 * nb.range

for(i in seq_along(nb.range)) {

nbSteps <- nb.range[i]

dtSim <- seq(as.timeDate(dtStart), as.timeDate(dtExpiry),
             length.out=nbSteps+1)

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

expiry.error <- tail(res$wealth,1)

mean.error[i] <- mean(expiry.error)
sd.error[i] <- sd(as.vector(expiry.error)) 
}

x11()
plotCI(x=nb.range, y=mean.error, uiw=sd.error/2, liw=sd.error/2, ylim=c(min(mean.error-sd.error/2), max(mean.error+sd.error/2)), 
       main='Delta hedging error vs. hedging frequency', xlab='Nb of rebalancings', ylab='Hedging error')


q <- lm((mean.error+sd.error/2) ~ sqrt(1/nb.range))
lines(nb.range, q$fitted.values, lwd=2, col='red')

##
## Transaction cost, round trip
##

tr.range <- seq(0, 1, .1)/100
nbSteps <- 250

dtSim <- seq(as.timeDate(dtStart), as.timeDate(dtExpiry),
             length.out=nbSteps+1)

dt <- as.real(dtExpiry-dtStart)/(365*nbSteps)

# price paths 
tSpot <- pathSimulator(dtSim = dtSim, nbPaths=nbPaths, 
    innovations.gen=sobolInnovations, path.gen=logNormal, 
    path.param = list(mu=0, sigma=sigma), S0=100, antithetic = F, 
    standardization = TRUE, trace = F)

# derived environment for scenario analysis
sce.env <- DataProvider(parent=base.env)
setData(sce.env, underlying, 'Price',
        time(tSpot), as.matrix(tSpot))

mean.error <- 0.0 * tr.range
sd.error <- 0.0 * tr.range

for(i in seq_along(tr.range)) {

# simulate a delta-hedge strategy along each path

assets = list(a)
res <- deltaHedge(assets, sce.env,
                  params=list(dtSim=time(tSpot),
                  transaction.cost=tr.range[i]/2), trace=F)

expiry.error <- tail(res$wealth,1)

mean.error[i] <- mean(expiry.error)
sd.error[i] <- sd(as.vector(expiry.error)) 
}

x11()

##
## compute adjustment to vol by Leland's model 
## tr.cost is round trip: (ask-bid)/mid

sigma.leland <- sigma * sqrt(1+ sqrt(2/pi) * (tr.range/(sigma*sqrt(dt))))

sigma.sim <- sigma - mean.error/vega
 
# create extra margin room on the right for an axis 
par(mar=c(5, 4, 4, 8) + 0.1)

# plot x vs. y 
plotCI(x=tr.range, y=mean.error, uiw=sd.error/2, liw=sd.error/2, ylim=c(min(mean.error-sd.error/2), max(mean.error+sd.error/2)), 
       xlab='', ylab='')

par(new=T)
# add x vs. 1/x 
plot(x=tr.range, y=sigma.leland, type="b", axes=F, pch=22, col="blue", lty=2, xlab="", ylab="")
lines(x=tr.range, y=sigma.sim, type="b", axes=F, pch=23, col="green", lty=2, xlab="", ylab="")
axis(side=4, las=2, cex.axis=0.7, tck=-.01)

legend('top', c("Leland", "Simulation"), lty=c(1,1), col=c("blue", "green"))

# add a title for the right axis 
mtext("adj. vol", side=4, line=3, cex.lab=1,las=2)

# add a main title and bottom and left axis labels 
title("Hedging error vs. transaction cost", xlab="% round trip cost",
   ylab="Hedging error")


