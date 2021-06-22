
library(latentcor)

test_that("zratio is NULL if type is continuous", {
  X = GenData(types = "con")
  expect_true(is.null(unlist(zratios(X = X, types_code = 0))))
})

test_that("zratio is the proportion of zeros if type is binary", {
  X = GenData(types = "bin")
  expect_equal(unlist(zratios(X = X, types_code = 1)), colMeans(X == 0))
})

test_that("zratio is the proportion of zeros if type is truncated", {
  X = GenData(type = "tru")
  expect_equal(unlist(zratios(X = X, types_code = 2)), colMeans(X == 0))
})

test_that("zratio is the proportion of zeros and ones if type is ternary", {
  X = GenData(type = "ter", copulas = c(NA, NA), c(.3, .5))
  expect_equal(unlist(zratios(X = X, types_code = 3)), c(colMeans(X2 == 0), 1 - colMeans(X2 == 2)))
})
