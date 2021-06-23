
library(latentcor)

K = runif(1, -.5, .5)
zratio1 = .5; zratio2 = .5

test_that("r should not exceed one.", {
  expect_lte(r_ml(K = K, zratio1 = zratio1, zratio2 = NULL, comb = "10"), 1)
  expect_lte(r_ml(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = "11"), 1)
  expect_lte(r_ml(K = K, zratio1 = zratio1, zratio2 = NULL, comb = "20"), 1)
  expect_lte(r_ml(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = "21"), 1)
  expect_lte(r_ml(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = "22"), 1)
})

zratio1 = matrix(c(.3, .7), ncol = 1)

test_that("r should not exceed one.", {
  expect_lte(r_ml(K = K, zratio1 = zratio1, zratio2 = NULL, comb = "30"), 1)
  expect_lte(r_ml(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = "31"), 1)
  expect_lte(r_ml(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = "32"), 1)
})

zratio2 = matrix(c(.2, .8), ncol = 1)
test_that("r should not exceed one.", {
  expect_lte(r_ml(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = "33"), 1)
})

K = sort(c(runif(1, -1, 0), runif(1, 0, 1)))
zratio1 = .5; zratio2 = .5

test_that("r is increasing in tau.", {
  expect_lte(r_ml(K = K[1], zratio1 = zratio1, zratio2 = NULL, comb = "10"),
             r_ml(K = K[2], zratio1 = zratio1, zratio2 = NULL, comb = "10"))
  expect_lte(r_ml(K = K[1], zratio1 = zratio1, zratio2 = zratio2, comb = "11"),
             r_ml(K = K[2], zratio1 = zratio1, zratio2 = zratio2, comb = "11"))
  expect_lte(r_ml(K = K[1], zratio1 = zratio1, zratio2 = NULL, comb = "20"),
             r_ml(K = K[2], zratio1 = zratio1, zratio2 = NULL, comb = "20"))
  expect_lte(r_ml(K = K[1], zratio1 = zratio1, zratio2 = zratio2, comb = "21"),
             r_ml(K = K[2], zratio1 = zratio1, zratio2 = zratio2, comb = "21"))
  expect_lte(r_ml(K = K[1], zratio1 = zratio1, zratio2 = zratio2, comb = "22"),
             r_ml(K = K[2], zratio1 = zratio1, zratio2 = zratio2, comb = "22"))
})

zratio1 = matrix(c(.3, .7), ncol = 1)

test_that("r is increasing in tau.", {
  expect_lte(r_ml(K = K[1], zratio1 = zratio1, zratio2 = NULL, comb = "30"),
             r_ml(K = K[2], zratio1 = zratio1, zratio2 = NULL, comb = "30"))
  expect_lte(r_ml(K = K[1], zratio1 = zratio1, zratio2 = zratio2, comb = "31"),
             r_ml(K = K[2], zratio1 = zratio1, zratio2 = zratio2, comb = "31"))
  expect_lte(r_ml(K = K[1], zratio1 = zratio1, zratio2 = zratio2, comb = "32"),
             r_ml(K = K[2], zratio1 = zratio1, zratio2 = zratio2, comb = "32"))
})

zratio2 = matrix(c(.2, .8), ncol = 1)
test_that("r is increasing in tau.", {
  expect_lte(r_ml(K = K[1], zratio1 = zratio1, zratio2 = zratio2, comb = "33"),
             r_ml(K = K[2], zratio1 = zratio1, zratio2 = zratio2, comb = "33"))
})
