
library(latentcor)

tau = runif(1, -1, 1)
zratio1 = .5; zratio2 = .5

test_that("r should not exceed one.", {
  expect_lte(r_sol(type1 = "binary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, tol = 10^(-6)), 1)
  expect_lte(r_sol(type1 = "binary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 10^(-6)), 1)
  expect_lte(r_sol(type1 = "trunc", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, tol = 10^(-6)), 1)
  expect_lte(r_sol(type1 = "trunc", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 10^(-6)), 1)
  expect_lte(r_sol(type1 = "trunc", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 10^(-6)), 1)
})

zratio1 = matrix(c(.3, .7), nrow = 1)

test_that("r should not exceed one.", {
  expect_lte(r_sol(type1 = "ternary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, tol = 10^(-6)), 1)
  expect_lte(r_sol(type1 = "ternary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 10^(-6)), 1)
  expect_lte(r_sol(type1 = "ternary", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 10^(-6)), 1)
})

zratio2 = matrix(c(.2, .8), nrow = 1)
test_that("r should not exceed one.", {
  expect_lte(r_sol(type1 = "ternary", type2 = "ternary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 10^(-6)), 1)
})
