
library(latentcor)

### Data setting
n = 100; p = 2 # sample size and dimensions for two datasets.

test_that("autocor is symmetric.", {
  expect_equal(autocor(p, 0.7), t(autocor(p, 0.7)))
})

Sigma = autocor(p, 0.7)
mu = rbinom(p, 1, 0.5)
dat = MASS::mvrnorm(n, mu = mu, Sigma = Sigma) # generate a data matrix of size.

test_that("truncated data have lower bound zero.", {
  expect_equal(apply(fromZtoX(Z = dat, copula = "cube", type = "trunc", c = matrix(rep(0, p), ncol = p), q = NULL), 2, min), rep(0, p))
})

test_that("binary data either zero or one.", {
  expect_equal(sort(unique(c(fromZtoX(Z = dat, copula = "cube", type = "binary", c = matrix(rep(0, p), ncol = p), q = NULL)))), c(0, 1))
})

test_that("ternary data should be zero, one or two.", {
  expect_equal(sort(unique(c(fromZtoX(Z = dat, copula = "cube", type = "ternary", c = matrix(rep(0:1, p), ncol = p), q = NULL)))), c(0, 1, 2))
})

p1 = 1; p2 = 2
Sigma = autocor(p1 + p2, 0.7)
mu = rbinom(p1 + p2, 1, 0.5)

test_that("truncated data have lower bound zero.", {
  expect_equal(apply(GenData(n = n, copula1 = "cube", copula2 = "cube", type1 = "ternary", type2 = "trunc", mu = mu, Sigma = Sigma, p1 = p1, p2 = p2,c1 = matrix(rep(0:1, p1), nrow = 2, ncol = p1), c2 = matrix(rep(0, p2), ncol = p2))$X2, 2, min), rep(0, p2))
})

test_that("truncated data have lower bound zero.", {
  expect_equal(sort(unique(c(GenData(n = n, copula1 = "cube", copula2 = "cube", type1 = "ternary", type2 = "binary", mu = mu, Sigma = Sigma, p1 = p1, p2 = p2,c1 = matrix(rep(0:1, p1), nrow = 2, ncol = p1), c2 = matrix(rep(0, p2), ncol = p2))$X2))), c(0, 1))
})

test_that("ternary data should be zero, one or two.", {
  expect_equal(sort(unique(c(GenData(n = n, copula1 = "cube", copula2 = "cube", type1 = "ternary", type2 = "trunc", mu = mu, Sigma = Sigma, p1 = p1, p2 = p2,c1 = matrix(rep(0:1, p1), nrow = 2, ncol = p1), c2 = matrix(rep(0, p2), ncol = p2))$X1))), c(0, 1, 2))
})


