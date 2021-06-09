
library(latentcor)

zratio1 = .5; zratio2 = .5

test_that("boundary for tau does not exceed one", {
  expect_lte(bound_bc(zratio1 = zratio1), 1)
  expect_lte(bound_bb(zratio1 = zratio1, zratio2 = zratio2), 1)
  expect_lte(bound_tc(zratio1 = zratio1), 1)
  expect_lte(bound_tb(zratio1 = zratio1, zratio2 = zratio2), 1)
  expect_lte(bound_tt(zratio1 = zratio1, zratio2 = zratio2), 1)
})

zratio1 = matrix(c(.3, .7), nrow = 1)
test_that("boundary for tau does not exceed one", {
  expect_lte(bound_nc(zratio1 = zratio1), 1)
  expect_lte(bound_nb(zratio1 = zratio1, zratio2 = zratio2), 1)
  expect_lte(bound_nt(zratio1 = zratio1, zratio2 = zratio2), 1)
})

zratio2 = matrix(c(.2, .8), nrow = 1)
test_that("boundary for tau does not exceed one", {
  expect_lte(bound_nn(zratio1 = zratio1, zratio2 = zratio2), 1)
})
