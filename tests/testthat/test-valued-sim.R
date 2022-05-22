library(ergm.count)
load("testnet3d.RData")

test_that("Poisson-reference ERGM with mean 2, TNT proposal", {
  set.seed(0)
  s <- simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", coef=log(2), output="stats", control=control.simulate.formula(MCMC.prop.args=list(p0=0.8)))
  test <- approx.hotelling.diff.test(s/6,mu0=2)
  expect_gt(test$p.value, 0.001)
})

test_that("Poisson-reference ERGM with mean 2, random proposal", {
  set.seed(0)
  s <- simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", coef=log(2), output="stats", control=control.simulate.formula(MCMC.prop.weights="random"))
  test <- approx.hotelling.diff.test(s/6,mu0=2)
  expect_gt(test$p.value, 0.001)
})

test_that("Poisson-reference ERGM with mean 2, zero-modified proposal", {
  set.seed(0)
  s <- simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", coef=log(2), output="stats", control=control.simulate.formula(MCMC.prop.weights="0inflated", MCMC.prop.args=list(p0=0.8)))
  test <- approx.hotelling.diff.test(s/6,mu0=2)
  expect_gt(test$p.value, 0.001)
})

## Geometric-reference
test_that("Geometric-reference ERGM with mean 2", {
  set.seed(0)
  s <- simulate(testnet3d~sum, nsim=1000, reference=~Geometric, response="w", coef=log(2/3), output="stats")
  test <- approx.hotelling.diff.test(s/6,mu0=2)
  expect_gt(test$p.value, 0.001)
})

## Binomial-reference
test_that("Binomial-reference ERGM with mean 5 trials and probability of success 0.4 for a mean 2", {
  set.seed(0)
  s <- simulate(testnet3d~sum, nsim=1000, reference=~Binomial(5), response="w", coef=log(.4/.6), output="stats")
  test <- approx.hotelling.diff.test(s/6,mu0=2)
  expect_gt(test$p.value, 0.001)
})

test_that("Poisson-reference ERGM with mean 2 and constraints, TNT proposal", {
  set.seed(0)
  s <- simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", constraints=~egocentric(I(c(FALSE,FALSE,TRUE))), coef=log(2), output="stats", control=control.simulate.formula(MCMC.prop.args=list(p0=0.5)))
  test <- approx.hotelling.diff.test(s/2,mu0=2)
  expect_gt(test$p.value, 0.001)
})