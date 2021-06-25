
library(latentcor)

K = runif(1, -1, 1)
zratio1 = .5; zratio2 = .5

test_that("r should not exceed one.", {
  expect_lte(r_sol(K = K, zratio1 = zratio1, zratio2 = NA, comb = "10", tol = 1e-8, ratio = .9), 1)
  expect_lte(r_sol(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = "11", tol = 1e-8, ratio = .9), 1)
  expect_lte(r_sol(K = K, zratio1 = zratio1, zratio2 = NA, comb = "20", tol = 1e-8, ratio = .9), 1)
  expect_lte(r_sol(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = "21", tol = 1e-8, ratio = .9), 1)
  expect_lte(r_sol(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = "22", tol = 1e-8, ratio = .9), 1)
})

zratio1 = matrix(c(.3, .7), ncol = 1)

test_that("r should not exceed one.", {
  expect_lte(r_sol(K = K, zratio1 = zratio1, zratio2 = NA, comb = "30", tol = 1e-8, ratio = .9), 1)
  expect_lte(r_sol(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = "31", tol = 1e-8, ratio = .9), 1)
  expect_lte(r_sol(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = "32", tol = 1e-8, ratio = .9), 1)
})

zratio2 = matrix(c(.2, .8), ncol = 1)
test_that("r should not exceed one.", {
  expect_lte(r_sol(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = "33", tol = 1e-8, ratio = .9), 1)
})

K = sort(c(runif(1, -1, 0), runif(1, 0, 1)))
zratio1 = .5; zratio2 = .5

test_that("r is increasing in tau.", {
  expect_lte(r_sol(K = K[1], zratio1 = zratio1, zratio2 = NA, comb = "10", tol = 1e-8, ratio = .9),
             r_sol(K = K[2], zratio1 = zratio1, zratio2 = NA, comb = "10", tol = 1e-8, ratio = .9))
  expect_lte(r_sol(K = K[1], zratio1 = zratio1, zratio2 = zratio2, comb = "11", tol = 1e-8, ratio = .9),
             r_sol(K = K[2], zratio1 = zratio1, zratio2 = zratio2, comb = "11", tol = 1e-8, ratio = .9))
  expect_lte(r_sol(K = K[1], zratio1 = zratio1, zratio2 = NA, comb = "20", tol = 1e-8, ratio = .9),
             r_sol(K = K[2], zratio1 = zratio1, zratio2 = NA, comb = "20", tol = 1e-8, ratio = .9))
  expect_lte(r_sol(K = K[1], zratio1 = zratio1, zratio2 = zratio2, comb = "21", tol = 1e-8, ratio = .9),
             r_sol(K = K[2], zratio1 = zratio1, zratio2 = zratio2, comb = "21", tol = 1e-8, ratio = .9))
  expect_lte(r_sol(K = K[1], zratio1 = zratio1, zratio2 = zratio2, comb = "22", tol = 1e-8, ratio = .9),
             r_sol(K = K[2], zratio1 = zratio1, zratio2 = zratio2, comb = "22", tol = 1e-8, ratio = .9))
})

zratio1 = matrix(c(.3, .7), ncol = 1)

test_that("r is increasing in tau.", {
  expect_lte(r_sol(K = K[1], zratio1 = zratio1, zratio2 = NA, comb = "30", tol = 1e-8, ratio = .9),
             r_sol(K = K[2], zratio1 = zratio1, zratio2 = NA, comb = "30", tol = 1e-8, ratio = .9))
  expect_lte(r_sol(K = K[1], zratio1 = zratio1, zratio2 = zratio2, comb = "31", tol = 1e-8, ratio = .9),
             r_sol(K = K[2], zratio1 = zratio1, zratio2 = zratio2, comb = "31", tol = 1e-8, ratio = .9))
  expect_lte(r_sol(K = K[1], zratio1 = zratio1, zratio2 = zratio2, comb = "32", tol = 1e-8, ratio = .9),
             r_sol(K = K[2], zratio1 = zratio1, zratio2 = zratio2, comb = "32", tol = 1e-8, ratio = .9))
})

zratio2 = matrix(c(.2, .8), ncol = 1)
test_that("r is increasing in tau.", {
  expect_lte(r_sol(K = K[1], zratio1 = zratio1, zratio2 = zratio2, comb = "33", tol = 1e-8, ratio = .9),
             r_sol(K = K[2], zratio1 = zratio1, zratio2 = zratio2, comb = "33", tol = 1e-8, ratio = .9))
})
