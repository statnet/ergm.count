#  File tests/testthat/test-valued-fit.R in package ergm.count, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2024 Statnet Commons
################################################################################
test_that("Poisson-reference ERGM fit", {
  set.seed(0)

  n <- 5

  m <- matrix(rpois(n^2,2),n,n)
  diag(m) <- 0
  y <- as.network(m, matrix.type="a", directed=TRUE, ignore.eval=FALSE, names.eval="w")

  truth <- log(sum(m)/n/(n-1))
  diag(m) <- NA

  efit <- ergm(y ~ sum, response="w", reference=~Poisson, verbose=TRUE, control=control.ergm(MCMLE.effectiveSize=128))

  summary(efit)
  true.llk <- sum(dpois(na.omit(c(m)), exp(coef(efit)), log=TRUE)) - sum(dpois(na.omit(c(m)), 1, log=TRUE))

  expect_lt(abs(coef(efit)-truth)/sqrt(vcov(efit, sources="estimation")), 3)
  expect_lt(abs(true.llk - logLik(efit))/sqrt(attr(logLik(efit),"vcov")), 3)
})
