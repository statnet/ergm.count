library(statnet.common)
opttest({
rm(list=ls())
library(ergm.count)
{
data(zach)

# Fit a binomial-reference ERGM.

zach.fit1 <- ergm(zach~nonzero+sum+nodefactor("role",base=2)+absdiffcat("faction.id"),
                  response="contexts", reference=~Binomial(8),
                  control=control.ergm(MCMLE.trustregion=1000))

mcmc.diagnostics(zach.fit1)

summary(zach.fit1)


}
},"zach.Rd")
