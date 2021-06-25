
library(latentcor)

Z = rnorm(100)

test_that("binary data either zero or one.", {
  expect_equal(sort(unique(c(fromZtoX(z = Z, type = "bin", copula = "cube", xp = .5)))), c(0, 1))
})

test_that("truncated data have lower bound zero.", {
  expect_equal(min(fromZtoX(z = Z, type = "tru", copula = "cube", xp = .5)), 0)
})

test_that("ternary data should be zero, one or two.", {
  expect_equal(sort(unique(c(fromZtoX(z = Z, type = "ter", copula = "cube", xp = c(.3, .5))))), c(0, 1, 2))
})

test_that("binary data have lower bound zero.", {
  expect_equal(sort(unique(c(GenData(types = "bin", copulas = "cube", XP = .5)))), c(0, 1))
})

test_that("truncated data have lower bound zero.", {
  expect_equal(min(GenData(types = "tru", copulas = "cube", XP = .5)), 0)
})

test_that("ternary data should be zero, one or two.", {
  expect_equal(sort(unique(c(GenData(types = "ter", copulas = "cube", XP = list(c(.3, .5)))))), c(0, 1, 2))
})
