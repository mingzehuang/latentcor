
library(latentcor)

test_that("n_x is 0 if no tie", {
  x = unique(rnorm(10))
  expect_equal(n_x(x = x, length(x)), 0)
})

test_that("n_x is not 1 if tie", {
  x = unique(rnorm(10)); x[3] = x[4]
  expect_equal(n_x(x = x, length(x)), 1)
})

# Data generation
test_that("latentcor is symmetric.", {
  X = gen_data(types = c("con", "con"), rhos = 1)$X
  expect_equal(latentcor(X = X, types = c("con", "con"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R,
               t(latentcor(X = X, types = c("con", "con"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R))
  expect_equal(latentcor(X = X, types = c("con", "con"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R,
               latentcor(X = X, method = "original", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R)
})

test_that("latentcor is symmetric.", {
  X = gen_data(types = c("con", "con"))$X; colnames(X) = c("var1", "var2")
  expect_equal(latentcor(X = X, types = c("con", "con"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R,
               t(latentcor(X = X, types = c("con", "con"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R))
  expect_equal(latentcor(X = X, types = c("con", "con"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R,
               latentcor(X = X, method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R)

})

test_that("latentcor is symmetric.", {
  X = gen_data(types = c("bin", "bin"), rhos = 1)$X
  expect_equal(latentcor(X = X, types = c("bin", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R,
               t(latentcor(X = X, types = c("bin", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R))
  expect_equal(latentcor(X = X, types = c("bin", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R,
               latentcor(X = X, method = "original", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R)
})

 test_that("latentcor is symmetric.", {
   X = gen_data(types = c("bin", "bin"))$X; colnames(X) = c("var1", "var2")
   expect_equal(latentcor(X = X, types = c("bin", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R,
                t(latentcor(X = X, types = c("bin", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R))
   expect_equal(latentcor(X = X, types = c("bin", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R,
                latentcor(X = X, method = "approx", nu = 0.5, tol = 1e-8, ratio = .9, showplot = FALSE)$R)
 })

test_that("latentcor is symmetric.", {
  X = gen_data(types = c("tru", "tru"), rhos = 1)$X
  expect_equal(latentcor(X = X, types = c("tru", "tru"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               t(latentcor(X = X, types = c("tru", "tru"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R))
  expect_equal(latentcor(X = X, types = c("tru", "tru"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               latentcor(X = X, method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R)
})

 test_that("latentcor is symmetric.", {
   X = gen_data(types = c("tru", "tru"))$X; colnames(X) = c("var1", "var2")
   expect_equal(latentcor(X = X, types = c("tru", "tru"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                t(latentcor(X = X, types = c("tru", "tru"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R))
   expect_equal(latentcor(X = X, types = c("tru", "tru"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                latentcor(X = X, method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R)
 })

test_that("latentcor is symmetric.", {
  X = gen_data(types = c("ter", "ter"), rhos = 1)$X
  expect_equal(latentcor(X = X, types = c("ter", "ter"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               t(latentcor(X = X, types = c("ter", "ter"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R))
  expect_equal(latentcor(X = X, types = c("ter", "ter"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               latentcor(X = X, method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R)
})

 test_that("latentcor is symmetric.", {
   X = gen_data(types = c("ter", "ter"))$X; colnames(X) = c("var1", "var2")
   expect_equal(latentcor(X = X, types = c("ter", "ter"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                t(latentcor(X = X, types = c("ter", "ter"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R))
   expect_equal(latentcor(X = X, types = c("ter", "ter"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                latentcor(X = X, method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R)
 })

test_that("latentcor is symmetric.", {
  X = gen_data(types = c("con", "bin"))$X
  expect_equal(latentcor(X = X, types = c("con", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               t(latentcor(X = X, types = c("con", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R))
  expect_equal(latentcor(X = X, types = c("con", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               latentcor(X = X, method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R)
})

 test_that("latentcor is symmetric.", {
   X = gen_data(types = c("con", "bin"))$X; colnames(X) = c("var1", "var2")
   expect_equal(latentcor(X = X, types = c("con", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                t(latentcor(X = X, types = c("con", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R))
   expect_equal(latentcor(X = X, types = c("con", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                latentcor(X = X, method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R)
 })

test_that("latentcor is symmetric.", {
  X = gen_data(types = c("tru", "bin"))$X
  expect_equal(latentcor(X = X, types = c("tru", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               t(latentcor(X = X, types = c("tru", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R))
  expect_equal(latentcor(X = X, types = c("tru", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               latentcor(X = X, method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R)
})

 test_that("latentcor is symmetric.", {
   X = gen_data(types = c("tru", "bin"))$X; colnames(X) = c("var1", "var2")
   expect_equal(latentcor(X = X, types = c("tru", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                t(latentcor(X = X, types = c("tru", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R))
   expect_equal(latentcor(X = X, types = c("tru", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                latentcor(X = X, method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R)

 })


test_that("latentcor is symmetric.", {
  X = gen_data(types = c("tru", "con"))$X
  expect_equal(latentcor(X = X, types = c("tru", "con"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               t(latentcor(X = X, types = c("tru", "con"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R))
  expect_equal(latentcor(X = X, types = c("tru", "con"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               latentcor(X = X, method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R)
})

 test_that("latentcor is symmetric.", {
   X = gen_data(types = c("tru", "con"))$X; colnames(X) = c("var1", "var2")
   expect_equal(latentcor(X = X, types = c("tru", "con"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                t(latentcor(X = X, types = c("tru", "con"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R))
   expect_equal(latentcor(X = X, types = c("tru", "con"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                latentcor(X = X, method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R)
 })

test_that("latentcor is symmetric.", {
  X = gen_data(types = c("ter", "con"))$X
  expect_equal(latentcor(X = X, types = c("ter", "con"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               t(latentcor(X = X, types = c("ter", "con"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R))
  expect_equal(latentcor(X = X, types = c("ter", "con"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               latentcor(X = X, method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R)
})

 test_that("latentcor is symmetric.", {
   X = gen_data(types = c("ter", "con"))$X; colnames(X) = c("var1", "var2")
   expect_equal(latentcor(X = X, types = c("ter", "con"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                t(latentcor(X = X, types = c("ter", "con"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R))
   expect_equal(latentcor(X = X, types = c("ter", "con"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                latentcor(X = X, method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R)
 })

test_that("latentcor is symmetric.", {
  X = gen_data(types = c("ter", "bin"))$X
  expect_equal(latentcor(X = X, types = c("ter", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               t(latentcor(X = X, types = c("ter", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R))
  expect_equal(latentcor(X = X, types = c("ter", "bin"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               latentcor(X = X, method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R)
})

 test_that("latentcor is symmetric.", {
   X = gen_data(types = c("ter", "bin"))$X; colnames(X) = c("var1", "var2")
   expect_equal(latentcor(X = X, types = c("ter", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                t(latentcor(X = X, types = c("ter", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R))
   expect_equal(latentcor(X = X, types = c("ter", "bin"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                latentcor(X = X, method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R)
 })

test_that("latentcor is symmetric.", {
  X = gen_data(types = c("ter", "tru"))$X
  expect_equal(latentcor(X = X, types = c("ter", "tru"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               t(latentcor(X = X, types = c("ter", "tru"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R))
  expect_equal(latentcor(X = X, types = c("ter", "tru"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               latentcor(X = X, method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R)
})

 test_that("latentcor is symmetric.", {
   X = gen_data(types = c("ter", "tru"))$X; colnames(X) = c("var1", "var2")
   expect_equal(latentcor(X = X, types = c("ter", "tru"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                t(latentcor(X = X, types = c("ter", "tru"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R))
   expect_equal(latentcor(X = X, types = c("ter", "tru"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
                latentcor(X = X, method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R)
 })

test_that("latentcor is symmetric.", {
  X = gen_data(types = c("ter", "ter"))$X
  expect_equal(latentcor(X = X, types = c("ter", "ter"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               t(latentcor(X = X, types = c("ter", "ter"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R))
  expect_equal(latentcor(X = X, types = c("ter", "ter"), method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               latentcor(X = X, method = "original", nu = 0.5, tol = 1e-8, ratio = .9)$R)
})

test_that("latentcor is symmetric.", {
  X = gen_data(types = c("ter", "ter"))$X; colnames(X) = c("var1", "var2")
  expect_equal(latentcor(X = X, types = c("ter", "ter"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               t(latentcor(X = X, types = c("ter", "ter"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R))
  expect_equal(latentcor(X = X, types = c("ter", "ter"), method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R,
               latentcor(X = X, method = "approx", nu = 0.5, tol = 1e-8, ratio = .9)$R)
})
