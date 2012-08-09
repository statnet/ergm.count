library(ergm.count)

## Poisson-reference
cat("Poisson-reference ERGM\n")
load("testnet3u.RData")

theta<-1
cat("Target mean:",exp(theta),"\n",sep="")

s<-simulate(testnet3u~sum, nsim=1000, reference="Poisson", response="w", coef=theta, statsonly=TRUE,
            control=control.simulate(MCMC.burnin=10000))

cat("Simulated mean (statsonly):",mean(s)/3,"\n",sep="")

s.full<-simulate(testnet3u~sum, nsim=1000, reference="Poisson", response="w", coef=theta, statsonly=FALSE,
            control=control.simulate(MCMC.burnin=10000))

cat("Simulated mean (full, computed):",mean(sapply(s.full,function(x) sum(x%e%"w")))/3,"\n",sep="")
cat("Simulated mean (full, stats):",mean(attr(s.full,"stats"))/3,"\n",sep="")

## Poisson-reference, zero-inflated
cat("Poisson-reference ERGM with zero-inflation\n")

theta<-1
cat("Target mean:",exp(theta),"\n",sep="")

s<-simulate(testnet3u~sum, nsim=1000, reference="Poisson", response="w", coef=theta, statsonly=TRUE,
            control=control.simulate(MCMC.burnin=10000, MCMC.prop.weights="0inflated"))

cat("Simulated mean (statsonly):",mean(s)/3,"\n",sep="")

s.full<-simulate(testnet3u~sum, nsim=1000, reference="Poisson", response="w", coef=theta, statsonly=FALSE,
                 control=control.simulate(MCMC.burnin=10000, MCMC.prop.weights="0inflated"))

cat("Simulated mean (full, computed):",mean(sapply(s.full,function(x) sum(x%e%"w")))/3,"\n",sep="")
cat("Simulated mean (full, stats):",mean(attr(s.full,"stats"))/3,"\n",sep="")
