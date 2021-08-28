
library(latentcor)

test_that("zratio is NULL if type is continuous", {
  X = gen_data(types = "con", copulas = "no")$X
  expect_true(is.na(unlist(zratios(X = X, types = "con"))))
})

test_that("zratio is the proportion of zeros if type is binary", {
  X = gen_data(types = "bin", copulas = "no")$X
  expect_equal(unlist(zratios(X = X, types = "bin")), as.numeric(colMeans(as.matrix(X == 0))))
})

test_that("zratio is the proportion of zeros if type is truncated", {
  X = gen_data(types = "tru", copulas = "no")$X
  expect_equal(unlist(zratios(X = X, types = "tru")), as.numeric(colMeans(as.matrix(X == 0))))
})

test_that("zratio is the proportion of zeros and ones if type is ternary", {
  X = gen_data(types = "ter", copulas = "no")$X
  expect_equal(unlist(zratios(X = X, types = "ter")), as.numeric(c(colMeans(as.matrix(X == 0)), 1 - colMeans(as.matrix(X == 2)))))
})
