#  File tests/testthat/test-valued-sim.R in package ergm.count, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2024 Statnet Commons
################################################################################
testnet3d <- network.initialize(3, directed=TRUE) # 6 dyads
testnet3u <- network.initialize(3, directed=FALSE) # 3 dyads

# Poisson-reference
test_that("Poisson-reference ERGM with mean 2, TNT proposal", {
  set.seed(0)
  ## expect_message(
    s <- simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", coef=log(2), output="stats", control=control.simulate.formula(MCMC.prop.args=list(p0=0.8)), verbose=TRUE)## ,
  ##   "\\bergm.count:MH_DiscTNT\\b"
  ## )
  test <- approx.hotelling.diff.test(s/6, mu0=2)
  expect_gt(test$p.value, 0.001)
})

test_that("Poisson-reference ERGM with mean 2 and constraints, TNT proposal", {
  set.seed(0)
  ## expect_message(
    s <- simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", constraints=~egocentric(I(c(FALSE,FALSE,TRUE))), coef=log(2), output="stats", control=control.simulate.formula(MCMC.prop.args=list(p0=0.5)), verbose=TRUE)## ,
  ##   "\\bergm.count:MH_DiscTNT\\b"
  ## )
  test <- approx.hotelling.diff.test(s/2, mu0=2)
  expect_gt(test$p.value, 0.001)
})

test_that("Poisson-reference ERGM with mean 2, random proposal", {
  set.seed(0)
  ## expect_message(
    s <- simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", coef=log(2), output="stats", control=control.simulate.formula(MCMC.prop.weights="random"), verbose=TRUE)## ,
  ##   "\\bergm.count:MH_Disc\\b"
  ## )
  test <- approx.hotelling.diff.test(s/6, mu0=2)
  expect_gt(test$p.value, 0.001)
})

test_that("Poisson-reference ERGM with mean 2, zero-modified proposal", {
  set.seed(0)
  ## expect_message(
    s <- simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", coef=log(2), output="stats", control=control.simulate.formula(MCMC.prop.weights="0inflated", MCMC.prop.args=list(p0=0.8)), verbose=TRUE)## ,
  ##   "\\bergm.count:MH_ZIPoisson\\b"
  ## )
  test <- approx.hotelling.diff.test(s/6, mu0=2)
  expect_gt(test$p.value, 0.001)
})

## Geometric-reference
test_that("Geometric-reference ERGM with mean 2, TNT proposal", {
  set.seed(0)
  ## expect_message(
    s <- simulate(testnet3d~sum, nsim=1000, reference=~Geometric, response="w", coef=log(2/3), output="stats", verbose=TRUE)## ,
  ##   "\\bergm.count:MH_DiscTNT\\b"
  ## )
  test <- approx.hotelling.diff.test(s/6, mu0=2)
  expect_gt(test$p.value, 0.001)
})

## Binomial-reference
test_that("Binomial-reference ERGM with mean 5 trials and probability of success 0.4 for a mean 2, TNT proposal", {
  set.seed(0)
  ## expect_message(
    s <- simulate(testnet3d~sum, nsim=1000, reference=~Binomial(5), response="w", coef=log(.4/.6), output="stats", verbose=TRUE)## ,
  ##   "\\bergm.count:MH_DiscTNT\\b"
  ## )
  test <- approx.hotelling.diff.test(s/6, mu0=2)
  expect_gt(test$p.value, 0.001)
})

## DiscUnif-reference
test_that("DiscUnif-reference ERGM with range [0,...,4] for a mean 2, TNT proposal", {
  set.seed(0)
  ## expect_message(
    s <- simulate(testnet3u~sum, nsim=1000, reference=~DiscUnif(0,4), response="w", coef=0, output="stats", verbose=TRUE)## ,
  ##   "\\bergm.count:MH_DiscTNT\\b"
  ## )
  expect_equal(min(s)/3, 0)
  expect_equal(max(s)/3, 4)
  test <- approx.hotelling.diff.test(s/3, mu0=2)
  expect_gt(test$p.value, 0.001)
})

test_that("DiscUnif-reference ERGM with range [1,...,3] for a mean 2, random proposal should be autodetected", {
  set.seed(0)
  ## expect_message(
    s <- simulate(testnet3u~sum, nsim=1000, reference=~DiscUnif(1,3), response="w", coef=0, output="stats", verbose=TRUE)## ,
  ##   "\\bergm.count:MH_Disc\\b"
  ## )
  expect_equal(min(s)/3, 1)
  expect_equal(max(s)/3, 3)
  test <- approx.hotelling.diff.test(s/3, mu0=2)
  expect_gt(test$p.value, 0.001)
})

test_that("DiscUnif-reference ERGM with range [-1,...,2] for a mean 1/2, TNT proposal", {
  set.seed(0)
  ## expect_message(
    s <- simulate(testnet3u~sum, nsim=1000, reference=~DiscUnif(-1,2), response="w", coef=0, output="stats", verbose=TRUE)## ,
  ##   "\\bergm.count:MH_DiscTNT\\b"
  ## )
  expect_equal(min(s)/3, -1)
  expect_equal(max(s)/3, 2)
  test <- approx.hotelling.diff.test(s/3, mu0=1/2)
  expect_gt(test$p.value, 0.001)
})

EDU <- function(theta, a, b){
  x <- a:b
  p <- exp(x*theta)
  sum(p*x)/sum(p)
}

test_that("DiscUnif-reference ERGM with range [-1,...,1] and a truncated geometric, TNT proposal", {
  set.seed(0)
  ## expect_message(
    s <- simulate(testnet3u~sum, nsim=1000, reference=~DiscUnif(-1,1), response="w", coef=1, output="stats", verbose=TRUE)## ,
  ##   "\\bergm.count:MH_DiscTNT\\b"
  ## )
  expect_equal(min(s)/3, -1)
  expect_equal(max(s)/3, 1)
  test <- approx.hotelling.diff.test(s/3, mu0=EDU(1, -1, 1))
  expect_gt(test$p.value, 0.001)
})

test_that("DiscUnif-reference ERGM with range [0,...,2] and a truncated geometric, TNT proposal", {
  set.seed(0)
  ## expect_message(
    s <- simulate(testnet3u~sum, nsim=1000, reference=~DiscUnif(0,2), response="w", coef=1, output="stats", verbose=TRUE)## ,
  ##   "\\bergm.count:MH_DiscTNT\\b"
  ## )
  expect_equal(min(s)/3, 0)
  expect_equal(max(s)/3, 2)
  test <- approx.hotelling.diff.test(s/3, mu0=EDU(1, 0, 2))
  expect_gt(test$p.value, 0.001)
})

