
library(latentcor)

test_that("zratio is NULL if type is continuous", {
  X = GenData(types = "con", copulas = NA)
  expect_true(is.null(unlist(zratios(X = X, types_code = 0))))
})

test_that("zratio is the proportion of zeros if type is binary", {
  X = GenData(types = "bin", copulas = NA, XP = .5)
  expect_equal(unlist(zratios(X = X, types_code = 1)), colMeans(as.matrix(X == 0)))
})

test_that("zratio is the proportion of zeros if type is truncated", {
  X = GenData(types = "tru", copulas = NA, XP = .5)
  expect_equal(unlist(zratios(X = X, types_code = 2)), colMeans(as.matrix(X == 0)))
})

test_that("zratio is the proportion of zeros and ones if type is ternary", {
  X = GenData(types = "ter", copulas = NA, XP = c(.3, .5))
  expect_equal(unlist(zratios(X = X, types_code = 3)), c(colMeans(as.matrix(X == 0)), 1 - colMeans(as.matrix(X == 2))))
})
