
library(latentcor)

### Data setting
n = 1000; p1 = 3; p2 = 2 # sample size and dimensions for two datasets.
rho = .9
# Data generation
simdata = GenData(n = n, type1 = "continuous", type2 = "binary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = NULL, c2 = 0)
X1 = simdata$X1; X2 = simdata$X2

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "continuous", method = "original", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X1, type1 = "continuous", method = "original", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "continuous", method = "ml", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X1, type1 = "continuous", method = "ml", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X2, type1 = "binary", method = "original", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "binary", method = "original", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X2, type1 = "binary", method = "ml", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "binary", method = "ml", tol = 1e-8, ratio = .9)))
})

simdata = GenData(n = n, type1 = "trunc", type2 = "ternary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = 0, c2 =  c(0, 1))
X1 = simdata$X1; X2 = simdata$X2

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "trunc", method = "original", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X1, type1 = "trunc", method = "original", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "trunc", method = "ml", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X1, type1 = "trunc", method = "ml", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X2, type1 = "ternary", method = "original", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "ternary", method = "original", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X2, type1 = "ternary", method = "ml", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "ternary", method = "ml", tol = 1e-8, ratio = .9)))
})

simdata = GenData(n = n, type1 = "continuous", type2 = "binary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = NULL, c2 = 0)
X1 = simdata$X1; X2 = simdata$X2

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "continuous", X2 = X2, type2 = "binary", method = "original", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "binary", X2 = X1, type2 = "continuous", method = "original", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "continuous", X2 = X2, type2 = "binary", method = "ml", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "binary", X2 = X1, type2 = "continuous", method = "ml", tol = 1e-8, ratio = .9)))
})

simdata = GenData(n = n, type1 = "trunc", type2 = "binary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = 0, c2 = 0)
X1 = simdata$X1; X2 = simdata$X2

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "trunc", X2 = X2, type2 = "binary", method = "original", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "binary", X2 = X1, type2 = "trunc", method = "original", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "trunc", X2 = X2, type2 = "binary", method = "ml", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "binary", X2 = X1, type2 = "trunc", method = "ml", tol = 1e-8, ratio = .9)))
})

simdata = GenData(n = n, type1 = "trunc", type2 = "continuous", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = 0, c2 = NULL)
X1 = simdata$X1; X2 = simdata$X2

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "trunc", X2 = X2, type2 = "continuous", method = "original", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "continuous", X2 = X1, type2 = "trunc", method = "original", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "trunc", X2 = X2, type2 = "continuous", method = "ml", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "continuous", X2 = X1, type2 = "trunc", method = "ml", tol = 1e-8, ratio = .9)))
})

simdata = GenData(n = n, type1 = "ternary", type2 = "continuous", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = c(0, 1), c2 =  NULL)
X1 = simdata$X1; X2 = simdata$X2

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "continuous", method = "original", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "continuous", X2 = X1, type2 = "ternary", method = "original", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "continuous", method = "ml", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "continuous", X2 = X1, type2 = "ternary", method = "ml", tol = 1e-8, ratio = .9)))
})

simdata = GenData(n = n, type1 = "ternary", type2 = "binary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = c(0, 1), c2 = 0)
X1 = simdata$X1; X2 = simdata$X2

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "binary", method = "original", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "binary", X2 = X1, type2 = "ternary", method = "original", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "binary", method = "ml", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "binary", X2 = X1, type2 = "ternary", method = "ml", tol = 1e-8, ratio = .9)))
})

simdata = GenData(n = n, type1 = "ternary", type2 = "trunc", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = c(0, 1), c2 = 0)
X1 = simdata$X1; X2 = simdata$X2

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "trunc", method = "original", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "trunc", X2 = X1, type2 = "ternary", method = "original", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "trunc", method = "ml", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "trunc", X2 = X1, type2 = "ternary", method = "ml", tol = 1e-8, ratio = .9)))
})

simdata = GenData(n = n, type1 = "ternary", type2 = "ternary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = c(0, 1), c2 = c(0, 1))
X1 = simdata$X1; X2 = simdata$X2

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "ternary", method = "original", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "ternary", X2 = X1, type2 = "ternary", method = "original", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "ternary", method = "ml", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "ternary", X2 = X1, type2 = "ternary", method = "ml", tol = 1e-8, ratio = .9)))
})

simdata = GenData(n = n, type1 = "continuous", type2 = "continuous", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = NULL, c2 = NULL)
X1 = simdata$X1; X2 = simdata$X2

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "continuous", X2 = X2, type2 = "continuous", method = "original", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "continuous", X2 = X1, type2 = "continuous", method = "original", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "continuous", X2 = X2, type2 = "continuous", method = "ml", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "continuous", X2 = X1, type2 = "continuous", method = "ml", tol = 1e-8, ratio = .9)))
})

simdata = GenData(n = n, type1 = "binary", type2 = "binary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = 0, c2 = 0)
X1 = simdata$X1; X2 = simdata$X2

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "binary", X2 = X2, type2 = "binary", method = "original", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "binary", X2 = X1, type2 = "binary", method = "original", tol = 1e-8, ratio = .9)))
})

test_that("estimateR is symmetric.", {
  expect_equal(estimateR(X1 = X1, type1 = "binary", X2 = X2, type2 = "binary", method = "ml", tol = 1e-8, ratio = .9),
               t(estimateR(X1 = X2, type1 = "binary", X2 = X1, type2 = "binary", method = "ml", tol = 1e-8, ratio = .9)))
})

simdata = GenData(n = n, type1 = "trunc", type2 = "trunc", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = 0, c2 = 0)
X1 = simdata$X1; X2 = simdata$X2

test_that("estimateR is symmetric.", {
  expect_equal(round(estimateR(X1 = X1, type1 = "trunc", X2 = X2, type2 = "trunc", method = "original", tol = 1e-8, ratio = .9), 4),
               t(round(estimateR(X1 = X2, type1 = "trunc", X2 = X1, type2 = "trunc", method = "original", tol = 1e-8, ratio = .9), 4)))
})

test_that("estimateR is symmetric.", {
  expect_equal(round(estimateR(X1 = X1, type1 = "trunc", X2 = X2, type2 = "trunc", method = "ml", tol = 1e-8, ratio = .9), 4),
               t(round(estimateR(X1 = X2, type1 = "trunc", X2 = X1, type2 = "trunc", method = "ml", tol = 1e-8, ratio = .9), 4)))
})
