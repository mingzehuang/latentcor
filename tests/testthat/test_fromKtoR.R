
library(latentcor)

p1 = 1; p2 = 1
Sigma = autocor(p1 + p2, 0.3)

zratio1 = as.matrix(c(.5, .6), ncol = 1)

test_that("fromKtoR should have consistent methods inside.", {
  expect_equal(fromKtoR(K = Sigma, zratio1 = NULL, type1 = "continuous", method = "original", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = NULL, type1 = "continuous", method = "approx", tol = 1e-8, ratio = -1))
})

test_that("fromKtoR should have consistent methods inside.", {
  expect_equal(fromKtoR(K = Sigma, zratio1 = NULL, type1 = "continuous", method = "ml", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = NULL, type1 = "continuous", method = "approx", tol = 1e-8, ratio = 2))
})

zratio1 = matrix(c(.3, .4, .6, .7), nrow = 2, ncol = 2)

test_that("fromKtoR should have consistent methods inside.", {
expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, type1 = "ternary", method = "original", tol = 1e-8, ratio = .9),
             fromKtoR(K = Sigma, zratio1 = zratio1, type1 = "ternary", method = "approx", tol = 1e-8, ratio = -1))
})

test_that("fromKtoR should have consistent methods inside.", {
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, type1 = "ternary", method = "ml", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, type1 = "ternary", method = "approx", tol = 1e-8, ratio = 2))
})


Sigma = matrix(c(.1, .2, .3, .4), 2, 2)
zratio1 = as.matrix(c(.5, .6), ncol = 1); zratio2 = as.matrix(c(.4, .7), ncol = 1)

test_that("fromKtoR should have consistent methods inside.", {
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "binary", type2 = "continuous", method = "original", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "binary", type2 = "continuous", method = "approx", tol = 1e-8, ratio = - 1))
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "binary", type2 = "binary", method = "original", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "binary", type2 = "binary", method = "approx", tol = 1e-8, ratio = - 1))
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "continuous", method = "original", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "continuous", method = "approx", tol = 1e-8, ratio = - 1))
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "binary", method = "original", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "binary", method = "approx", tol = 1e-8, ratio = - 1))
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "trunc", method = "original", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "trunc", method = "approx", tol = 1e-8, ratio = - 1))
})

test_that("fromKtoR should have consistent methods inside.", {
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "binary", type2 = "continuous", method = "ml", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "binary", type2 = "continuous", method = "approx", tol = 1e-8, ratio = 2))
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "binary", type2 = "binary", method = "ml", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "binary", type2 = "binary", method = "approx", tol = 1e-8, ratio = 2))
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "continuous", method = "ml", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "continuous", method = "approx", tol = 1e-8, ratio = 2))
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "binary", method = "ml", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "binary", method = "approx", tol = 1e-8, ratio = 2))
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "trunc", method = "ml", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "trunc", method = "approx", tol = 1e-8, ratio = 2))
})


zratio1 = matrix(c(.3, .4, .6, .7), nrow = 2, ncol = 2)

test_that("fromKtoR should have consistent methods inside.", {
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "continuous", method = "original", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "continuous", method = "approx", tol = 1e-8, ratio = -1))
})

test_that("fromKtoR should have consistent methods inside.", {
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "continuous", method = "ml", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "continuous", method = "approx", tol = 1e-8, ratio = 2))
})

test_that("fromKtoR should have consistent methods inside.", {
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "binary", method = "original", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "binary", method = "approx", tol = 1e-8, ratio = -1))
})

test_that("fromKtoR should have consistent methods inside.", {
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "binary", method = "ml", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "binary", method = "approx", tol = 1e-8, ratio = 2))
})

test_that("fromKtoR should have consistent methods inside.", {
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "trunc", method = "original", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "trunc", method = "approx", tol = 1e-8, ratio = -1))
})

test_that("fromKtoR should have consistent methods inside.", {
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "trunc", method = "ml", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "trunc", method = "approx", tol = 1e-8, ratio = 2))
})

zratio2 = matrix(c(.2, .3, .5, .8), nrow = 2, ncol = 2)

test_that("fromKtoR should have consistent methods inside.", {
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "ternary", method = "original", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "ternary", method = "approx", tol = 1e-8, ratio = -1))
})

test_that("fromKtoR should have consistent methods inside.", {
  expect_equal(fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "ternary", method = "ml", tol = 1e-8, ratio = .9),
               fromKtoR(K = Sigma, zratio1 = zratio1, zratio2 = zratio2, type1 = "ternary", type2 = "ternary", method = "approx", tol = 1e-8, ratio = 2))
})
