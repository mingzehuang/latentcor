
library(latentcor)

zratio1 = .5; zratio2 = .5

test_that("cutoff selection should go to original method if tau equals one.", {
  expect_true(cutoff(type1 = NULL, type2 = NULL, tau = 1, zratio1 = NULL, zratio2 = NULL, method = "original", ratio = NULL))
  expect_false(cutoff(type1 = NULL, type2 = NULL, tau = 1, zratio1 = NULL, zratio2 = NULL, method = "ml", ratio = NULL))
  expect_true(cutoff(type1 = "binary", type2 = "continuous", tau = 1, zratio1 = zratio1, zratio2 = NULL, method = "approx", ratio = .9))
  expect_true(cutoff(type1 = "binary", type2 = "binary", tau = 1, zratio1 = zratio1, zratio2 = zratio2, method = "approx", ratio = .9))
  expect_true(cutoff(type1 = "trunc", type2 = "continuous", tau = 1, zratio1 = zratio1, zratio2 = NULL, method = "approx", ratio = .9))
  expect_true(cutoff(type1 = "trunc", type2 = "binary", tau = 1, zratio1 = zratio1, zratio2 = zratio2, method = "approx", ratio = .9))
  expect_true(cutoff(type1 = "trunc", type2 = "trunc", tau = 1, zratio1 = zratio1, zratio2 = zratio2, method = "approx", ratio = .9))
})

zratio1 = matrix(c(.3, .7), nrow = 1)

test_that("cutoff selection should go to original method if tau equals one.", {
  expect_true(cutoff(type1 = "ternary", type2 = "continuous", tau = 1, zratio1 = zratio1, zratio2 = NULL, method = "approx", ratio = .9))
  expect_true(cutoff(type1 = "ternary", type2 = "binary", tau = 1, zratio1 = zratio1, zratio2 = zratio2, method = "approx", ratio = .9))
  expect_true(cutoff(type1 = "ternary", type2 = "trunc", tau = 1, zratio1 = zratio1, zratio2 = zratio2, method = "approx", ratio = .9))
})

zratio2 = matrix(c(.2, .8), nrow = 1)

test_that("cutoff selection should go to original method if tau equals one.", {
  expect_true(cutoff(type1 = "ternary", type2 = "ternary", tau = 1, zratio1 = zratio1, zratio2 = zratio2, method = "approx", ratio = .9))
})
