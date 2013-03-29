library(ergm.count)
opttest({
load("testnet3d.RData")

## Poisson-reference
cat("======== Poisson-reference ERGM with mean 2\n")
cat("==== statsonly=TRUE\n")
s <- simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", coef=log(2), statsonly=TRUE)
test <- approx.hotelling.diff.test(s/6,mu0=2)
if(test$p.value<.01) {print(test); stop("Simulation test failed.")}

cat("==== statsonly=FALSE\n")
s.full<-simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", coef=log(2), statsonly=FALSE)
test <- approx.hotelling.diff.test(sapply(s.full, function(x) sum(as.matrix(x,m="a",a="w"))/6,simplify=TRUE),mu0=2)
if(test$p.value<.01) {print(test); stop("Simulation test failed.")}

cat("======== Poisson-reference ERGM with mean 2, zero-modified proposal\n")
cat("==== statsonly=TRUE\n")
s <- simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", coef=log(2), statsonly=TRUE, control=control.simulate.formula(MCMC.prop.weights="0inflated", MCMC.prop.args=list(p0=0.8)))
test <- approx.hotelling.diff.test(s/6,mu0=2)
if(test$p.value<.01) {print(test); stop("Simulation test failed.")}

cat("==== statsonly=FALSE\n")
s.full<-simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", coef=log(2), statsonly=FALSE, control=control.simulate.formula(MCMC.prop.weights="0inflated", MCMC.prop.args=list(p0=0.8)))
test <- approx.hotelling.diff.test(sapply(s.full, function(x) sum(as.matrix(x,m="a",a="w"))/6,simplify=TRUE),mu0=2)
if(test$p.value<.01) {print(test); stop("Simulation test failed.")}

## Geometric-reference
cat("======== Geometric-reference ERGM with mean 2\n")
cat("==== statsonly=TRUE\n")
s <- simulate(testnet3d~sum, nsim=1000, reference=~Geometric, response="w", coef=log(2/3), statsonly=TRUE)
test <- approx.hotelling.diff.test(s/6,mu0=2)
if(test$p.value<.01) {print(test); stop("Simulation test failed.")}

cat("==== statsonly=FALSE\n")
s.full<-simulate(testnet3d~sum, nsim=1000, reference=~Geometric, response="w", coef=log(2/3), statsonly=FALSE)
test <- approx.hotelling.diff.test(sapply(s.full, function(x) sum(as.matrix(x,m="a",a="w"))/6,simplify=TRUE),mu0=2)
if(test$p.value<.01) {print(test); stop("Simulation test failed.")}

## Binomial-reference
cat("======== Binomial-reference ERGM with mean 5 trials and probability of success 0.4 for a mean 2\n")
cat("==== statsonly=TRUE\n")
s <- simulate(testnet3d~sum, nsim=1000, reference=~Binomial(5), response="w", coef=log(.4/.6), statsonly=TRUE)
test <- approx.hotelling.diff.test(s/6,mu0=2)
if(test$p.value<.01) {print(test); stop("Simulation test failed.")}

cat("==== statsonly=FALSE\n")
s.full<-simulate(testnet3d~sum, nsim=1000, reference=~Binomial(5), response="w", coef=log(.4/.6), statsonly=FALSE)
test <- approx.hotelling.diff.test(sapply(s.full, function(x) sum(as.matrix(x,m="a",a="w"))/6,simplify=TRUE),mu0=2)
if(test$p.value<.01) {print(test); stop("Simulation test failed.")}
})
