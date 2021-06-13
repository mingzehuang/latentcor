
library(latentcor)

### Data setting
n = 1000; p1 = 3; p2 = 2 # sample size and dimensions for two datasets.
rho = .9
# Data generation
simdata = GenData(n = n, type1 = "continuous", type2 = "binary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = NULL, c2 = 0)
X1 = simdata$X1; X2 = simdata$X2
colnames(X1) = c("X11", "X12", "X13"); colnames(X2) = c("X21", "X22")

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "continuous", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9),
               t(estR(X1 = X1, type1 = "continuous", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "continuous", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9),
               t(estR(X1 = X1, type1 = "continuous", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X2, type1 = "binary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9),
               t(estR(X1 = X2, type1 = "binary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X2, type1 = "binary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9),
               t(estR(X1 = X2, type1 = "binary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)))
})

simdata = GenData(n = n, type1 = "trunc", type2 = "ternary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = 0, c2 = c(0, 1))
X1 = simdata$X1; X2 = simdata$X2
colnames(X1) = c("X11", "X12", "X13"); colnames(X2) = c("X21", "X22")

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "trunc", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9),
               t(estR(X1 = X1, type1 = "trunc", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "trunc", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9),
               t(estR(X1 = X1, type1 = "trunc", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X2, type1 = "ternary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9),
               t(estR(X1 = X2, type1 = "ternary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X2, type1 = "ternary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9),
               t(estR(X1 = X2, type1 = "ternary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)))
})

simdata = GenData(n=n, type1 = "continuous", type2 = "binary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = NULL, c2 = 0)
X1 = simdata$X1; X2 = simdata$X2
colnames(X1) = c("X11", "X12", "X13"); colnames(X2) = c("X21", "X22")

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "continuous", X2 = X2, type2 = "binary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "binary", X2 = X1, type2 = "continuous", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "continuous", X2 = X2, type2 = "binary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "binary", X2 = X1, type2 = "continuous", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

simdata = GenData(n = n, type1 = "trunc", type2 = "binary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = 0, c2 = 0)
X1 = simdata$X1; X2 = simdata$X2
colnames(X1) = c("X11", "X12", "X13"); colnames(X2) = c("X21", "X22")

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "trunc", X2 = X2, type2 = "binary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "binary", X2 = X1, type2 = "trunc", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "trunc", X2 = X2, type2 = "binary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "binary", X2 = X1, type2 = "trunc", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

simdata = GenData(n = n, type1 = "trunc", type2 = "continuous", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = 0, c2 = NULL)
X1 = simdata$X1; X2 = simdata$X2
colnames(X1) = c("X11", "X12", "X13"); colnames(X2) = c("X21", "X22")

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "trunc", X2 = X2, type2 = "continuous", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "continuous", X2 = X1, type2 = "trunc", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "trunc", X2 = X2, type2 = "continuous", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "continuous", X2 = X1, type2 = "trunc", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

simdata = GenData(n = n, type1 = "ternary", type2 = "continuous", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = c(0, 1), c2 = NULL)
X1 = simdata$X1; X2 = simdata$X2
colnames(X1) = c("X11", "X12", "X13"); colnames(X2) = c("X21", "X22")

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "continuous", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "continuous", X2 = X1, type2 = "ternary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "continuous", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "continuous", X2 = X1, type2 = "ternary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

simdata = GenData(n = n, type1 = "ternary", type2 = "binary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = c(0, 1), c2 = 0)
X1 = simdata$X1; X2 = simdata$X2
colnames(X1) = c("X11", "X12", "X13"); colnames(X2) = c("X21", "X22")

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "binary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "binary", X2 = X1, type2 = "ternary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "binary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "binary", X2 = X1, type2 = "ternary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

simdata = GenData(n = n, type1 = "ternary", type2 = "trunc", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = c(0, 1), c2 = 0)
X1 = simdata$X1; X2 = simdata$X2
colnames(X1) = c("X11", "X12", "X13"); colnames(X2) = c("X21", "X22")

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "trunc", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "trunc", X2 = X1, type2 = "ternary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "trunc", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "trunc", X2 = X1, type2 = "ternary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

simdata = GenData(n = n, type1 = "ternary", type2 = "ternary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = c(0, 1), c2 = c(0, 1))
X1 = simdata$X1; X2 = simdata$X2
colnames(X1) = c("X11", "X12", "X13"); colnames(X2) = c("X21", "X22")

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "ternary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "ternary", X2 = X1, type2 = "ternary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "ternary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "ternary", X2 = X1, type2 = "ternary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

simdata = GenData(n = n, type1 = "continuous", type2 = "continuous", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = NULL, c2 = NULL)
X1 = simdata$X1; X2 = simdata$X2
colnames(X1) = c("X11", "X12", "X13"); colnames(X2) = c("X21", "X22")

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "continuous", X2 = X2, type2 = "continuous", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "continuous", X2 = X1, type2 = "continuous", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "continuous", X2 = X2, type2 = "continuous", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "continuous", X2 = X1, type2 = "continuous", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

simdata = GenData(n = n, type1 = "binary", type2 = "binary", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = 0, c2 = 0)
X1 = simdata$X1; X2 = simdata$X2
colnames(X1) = c("X11", "X12", "X13"); colnames(X2) = c("X21", "X22")

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "binary", X2 = X2, type2 = "binary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "binary", X2 = X1, type2 = "binary", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

test_that("estR is symmetric.", {
  expect_equal(estR(X1 = X1, type1 = "binary", X2 = X2, type2 = "binary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12,
               t(estR(X1 = X2, type1 = "binary", X2 = X1, type2 = "binary", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12))
})

simdata = GenData(n = n, type1 = "trunc", type2 = "trunc", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = 0, c2 = 0)
X1 = simdata$X1; X2 = simdata$X2
colnames(X1) = c("X11", "X12", "X13"); colnames(X2) = c("X21", "X22")

test_that("estR is symmetric.", {
  expect_equal(round(estR(X1 = X1, type1 = "trunc", X2 = X2, type2 = "trunc", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12, 2),
               t(round(estR(X1 = X2, type1 = "trunc", X2 = X1, type2 = "trunc", method = "original", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12, 2)))
})

test_that("estR is symmetric.", {
  expect_equal(round(estR(X1 = X1, type1 = "trunc", X2 = X2, type2 = "trunc", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12, 2),
               t(round(estR(X1 = X2, type1 = "trunc", X2 = X1, type2 = "trunc", method = "ml", use.nearPD = TRUE, nu = 0.5, tol = 1e-6, ratio = .9)$R12, 2)))
})
