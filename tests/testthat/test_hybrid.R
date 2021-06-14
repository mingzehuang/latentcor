
library(latentcor)
set.seed(0)
tau = runif(1, -1, 1)
zratio1 =as.matrix( .5); zratio2 =as.matrix( .5)

test_that("Method selection should be consistent with corresponding methods.", {
  expect_equal(R_sol(type1 = "binary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "orignal", tol = 1e-8, ratio = .9),
               r_sol(type1 = "binary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, tol = 1e-8))
  expect_equal(R_sol(type1 = "binary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "orignal", tol = 1e-8, ratio = .9),
               r_sol(type1 = "binary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 1e-8))
  expect_equal(R_sol(type1 = "trunc", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "orginal", tol = 1e-8, ratio = .9),
               r_sol(type1 = "trunc", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, tol = 1e-8))
  expect_equal(R_sol(type1 = "trunc", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "original", tol = 1e-8, ratio = .9),
               r_sol(type1 = "trunc", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 1e-8))
  expect_equal(R_sol(type1 = "trunc", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "original", tol = 1e-8, ratio = .9),
               r_sol(type1 = "trunc", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 1e-8))
})

test_that("Method selection should be consistent with corresponding methods.", {
  expect_equal(R_sol(type1 = "binary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "ml", tol = 1e-8, ratio = .9),
               r_ml(type1 = "binary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL))
  expect_equal(R_sol(type1 = "binary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "ml", tol = 1e-8, ratio = .9),
               r_ml(type1 = "binary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2))
  expect_equal(R_sol(type1 = "trunc", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "ml", tol = 1e-8, ratio = .9),
               r_ml(type1 = "trunc", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL))
  expect_equal(R_sol(type1 = "trunc", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "ml", tol = 1e-8, ratio = .9),
               r_ml(type1 = "trunc", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2))
  expect_equal(R_sol(type1 = "trunc", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "ml", tol = 1e-8, ratio = .9),
               r_ml(type1 = "trunc", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2))
})

test_that("Method selection should be consistent with corresponding methods.", {
  expect_equal(R_sol(type1 = "binary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "original", tol = 1e-8, ratio = .9),
               R_sol(type1 = "binary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "approx", tol = 1e-8, ratio = - 1))
  expect_equal(R_sol(type1 = "binary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "original", tol = 1e-8, ratio = .9),
               R_sol(type1 = "binary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "approx", tol = 1e-8, ratio = - 1))
  expect_equal(R_sol(type1 = "trunc", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "original", tol = 1e-8, ratio = .9),
               R_sol(type1 = "trunc", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "approx", tol = 1e-8, ratio = - 1))
  expect_equal(R_sol(type1 = "trunc", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "original", tol = 1e-8, ratio = .9),
               R_sol(type1 = "trunc", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "approx", tol = 1e-8, ratio = - 1))
  expect_equal(R_sol(type1 = "trunc", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "original", tol = 1e-8, ratio = .9),
               R_sol(type1 = "trunc", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "approx", tol = 1e-8, ratio = - 1))
})

test_that("Method selection should be consistent with corresponding methods.", {
  expect_equal(R_sol(type1 = "binary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "ml", tol = 1e-8, ratio = .9),
               R_sol(type1 = "binary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "approx", tol = 1e-8, ratio = 2))
  expect_equal(R_sol(type1 = "binary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "ml", tol = 1e-8, ratio = .9),
               R_sol(type1 = "binary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "approx", tol = 1e-8, ratio = 2))
  expect_equal(R_sol(type1 = "trunc", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "ml", tol = 1e-8, ratio = .9),
               R_sol(type1 = "trunc", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "approx", tol = 1e-8, ratio = 2))
  expect_equal(R_sol(type1 = "trunc", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "ml", tol = 1e-8, ratio = .9),
               R_sol(type1 = "trunc", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "approx", tol = 1e-8, ratio = 2))
  expect_equal(R_sol(type1 = "trunc", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "ml", tol = 1e-8, ratio = .9),
               R_sol(type1 = "trunc", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "approx", tol = 1e-8, ratio = 2))
})

zratio1 = matrix(c(.3, .7), nrow = 1)

test_that("Method selection should be consistent with corresponding methods.", {
  expect_equal(R_sol(type1 = "ternary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "original", tol = 1e-8, ratio = .9),
               r_sol(type1 = "ternary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, tol = 1e-8))
  expect_equal(R_sol(type1 = "ternary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "original", tol = 1e-8, ratio = .9),
               r_sol(type1 = "ternary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 1e-8))
  expect_equal(R_sol(type1 = "ternary", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "original", tol= 1e-8, ratio = .9),
               r_sol(type1 = "ternary", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 1e-8))
})

test_that("Method selection should be consistent with corresponding methods.", {
  expect_equal(R_sol(type1 = "ternary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "ml", tol = 1e-8, ratio = .9),
               r_ml(type1 = "ternary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL))
  expect_equal(R_sol(type1 = "ternary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "ml", tol = 1e-8, ratio = .9),
               r_ml(type1 = "ternary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2))
  expect_equal(R_sol(type1 = "ternary", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "ml", tol= 1e-8, ratio = .9),
               r_ml(type1 = "ternary", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2))
})

test_that("Method selection should be consistent with corresponding methods.", {
  expect_equal(R_sol(type1 = "ternary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "original", tol = 1e-8, ratio = .9),
               R_sol(type1 = "ternary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "approx", tol = 1e-8, ratio = - 1))
  expect_equal(R_sol(type1 = "ternary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "original", tol = 1e-8, ratio = .9),
               R_sol(type1 = "ternary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "approx", tol = 1e-8, ratio = - 1))
  expect_equal(R_sol(type1 = "ternary", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "original", tol= 1e-8, ratio = .9),
               R_sol(type1 = "ternary", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "approx", tol= 1e-8, ratio = - 1))
})

test_that("Method selection should be consistent with corresponding methods.", {
  expect_equal(R_sol(type1 = "ternary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "ml", tol = 1e-8, ratio = .9),
               R_sol(type1 = "ternary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, method = "approx", tol = 1e-8, ratio = 2))
  expect_equal(R_sol(type1 = "ternary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "ml", tol = 1e-8, ratio = .9),
               R_sol(type1 = "ternary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "approx", tol = 1e-8, ratio = 2))
  expect_equal(R_sol(type1 = "ternary", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "ml", tol= 1e-8, ratio = .9),
               R_sol(type1 = "ternary", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "approx", tol= 1e-8, ratio = 2))
})

zratio2 = matrix(c(.2, .8), nrow = 1)
test_that("Method selection should be consistent with corresponding methods.", {
  expect_equal(R_sol(type1 = "ternary", type2 = "ternary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "original", tol = 1e-8, ratio = .9),
               r_sol(type1 = "ternary", type2 = "ternary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 1e-8))
})

test_that("Method selection should be consistent with corresponding methods.", {
  expect_equal(R_sol(type1 = "ternary", type2 = "ternary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "original", tol = 1e-8, ratio = .9),
               R_sol(type1 = "ternary", type2 = "ternary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "approx", tol = 1e-8, ratio = - 1))
})

test_that("Method selection should be consistent with corresponding methods.", {
  expect_equal(R_sol(type1 = "ternary", type2 = "ternary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "ml", tol = 1e-8, ratio = .9),
               R_sol(type1 = "ternary", type2 = "ternary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, method = "approx", tol = 1e-8, ratio = 2))
})

