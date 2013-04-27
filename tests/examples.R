#  File tests/examples.R in package ergm.count, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
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
