
library(latentcor)

test_that("n_x is 0 if no tie", {
  x = unique(rnorm(10))
  expect_equal(n_x(x = x, length(x)), 0)
})

# Data generation
test_that("estR is symmetric.", {
  X = GenData(types = c("con", "con")); colnames(X) = c("con", "con")
  expect_equal(estR(X = X, types = c("con", "con"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("con", "con"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("con", "con")); colnames(X) = c("con", "con")
  expect_equal(estR(X = X, types = c("con", "con"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("con", "con"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("bin", "bin"), XP = list(.5, .5)); colnames(X) = c("bin", "bin")
  expect_equal(estR(X = X, types = c("bin", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("bin", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("bin", "bin"), XP = list(.5, .5)); colnames(X) = c("bin", "bin")
  expect_equal(estR(X = X, types = c("bin", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("bin", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("tru", "tru"), XP = list(.5, .5)); colnames(X) = c("tru", "tru")
  expect_equal(estR(X = X, types = c("tru", "tru"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("tru", "tru"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("tru", "tru"), XP = list(.5, .5)); colnames(X) = c("tru", "tru")
  expect_equal(estR(X = X, types = c("tru", "tru"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("tru", "tru"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("ter", "ter"), XP = list(c(.3, .5), c(.2, .6))); colnames(X) = c("ter", "ter")
  expect_equal(estR(X = X, types = c("ter", "ter"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("ter", "ter"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("ter", "ter"), XP = list(c(.3, .5), c(.2, .6))); colnames(X) = c("ter", "ter")
  expect_equal(estR(X = X, types = c("ter", "ter"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("ter", "ter"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("con", "bin"), XP = list(NA, .5)); colnames(X) = c("con", "bin")
  expect_equal(estR(X = X, types = c("con", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("con", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("con", "bin"), XP = list(NA, .5)); colnames(X) = c("con", "bin")
  expect_equal(estR(X = X, types = c("con", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("con", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("tru", "bin"), XP = list(.5, .5)); colnames(X) = c("tru", "bin")
  expect_equal(estR(X = X, types = c("tru", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("tru", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("tru", "bin"), XP = list(.5, .5)); colnames(X) = c("tru", "bin")
  expect_equal(estR(X = X, types = c("tru", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("tru", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})


test_that("estR is symmetric.", {
  X = GenData(types = c("tru", "con"), XP = list(.5, NA)); colnames(X) = c("tru", "con")
  expect_equal(estR(X = X, types = c("tru", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("tru", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("tru", "con"), XP = list(.5, NA)); colnames(X) = c("tru", "con")
  expect_equal(estR(X = X, types = c("tru", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("tru", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("ter", "con"), XP = list(c(.3, .6), NA)); colnames(X) = c("ter", "con")
  expect_equal(estR(X = X, types = c("ter", "con"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("ter", "con"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("ter", "con"), XP = list(c(.3, .6), NA)); colnames(X) = c("ter", "con")
  expect_equal(estR(X = X, types = c("ter", "con"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("ter", "con"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("ter", "bin"), XP = list(c(.3, .6), .5)); colnames(X) = c("ter", "bin")
  expect_equal(estR(X = X, types = c("ter", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("ter", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("ter", "bin"), XP = list(c(.3, .6), .5)); colnames(X) = c("ter", "bin")
  expect_equal(estR(X = X, types = c("ter", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("ter", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("ter", "tru"), XP = list(c(.3, .6), .5)); colnames(X) = c("ter", "tru")
  expect_equal(estR(X = X, types = c("ter", "tru"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("ter", "tru"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("ter", "tru"), XP = list(c(.3, .6), .5)); colnames(X) = c("ter", "tru")
  expect_equal(estR(X = X, types = c("ter", "tru"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("ter", "tru"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("ter", "ter"), XP = list(c(.3, .6), c(.2, .5))); colnames(X) = c("ter", "ter")
  expect_equal(estR(X = X, types = c("ter", "ter"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("ter", "ter"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})

test_that("estR is symmetric.", {
  X = GenData(types = c("ter", "ter"), XP = list(c(.3, .6), c(.2, .5))); colnames(X) = c("ter", "ter")
  expect_equal(estR(X = X, types = c("ter", "ter"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = TRUE)$R,
               t(estR(X = X, types = c("ter", "ter"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, corplot = FALSE)$R))
})
