
library(latentcor)

r = sort(runif(2, -1, 1)); zratio1 = .5; zratio2 = .5

test_that("bridgeF_bc is increasing in r", {
  expect_lte(bridgeF_bc(r = r[1], zratio1 = zratio1), bridgeF_bc(r = r[2], zratio1 = zratio1))
})

test_that("bridgeF_bb is increasing in r", {
  expect_lte(bridgeF_bb(r = r[1], zratio1 = zratio1, zratio2 = zratio2), bridgeF_bb(r = r[2], zratio1 = zratio1, zratio2 = zratio2))
})

test_that("bridgeF_tc is increasing in r", {
  expect_lte(bridgeF_tc(r = r[1], zratio1 = zratio1), bridgeF_tc(r = r[2], zratio1 = zratio1))
})

test_that("bridgeF_tb is increasing in r", {
  expect_lte(bridgeF_tb(r = r[1], zratio1 = zratio1, zratio2 = zratio2), bridgeF_tb(r = r[2], zratio1 = zratio1, zratio2 = zratio2))
})

test_that("bridgeF_tt is increasing in r", {
  expect_lte(bridgeF_tt(r = r[1], zratio1 = zratio1, zratio2 = zratio2), bridgeF_tt(r = r[2], zratio1 = zratio1, zratio2 = zratio2))
})

zratio1 = matrix(c(.3, .7), nrow = 1)

test_that("bridgeF_nc is increasing in r", {
  expect_lte(bridgeF_nc(r = r[1], zratio1 = zratio1), bridgeF_nc(r = r[2], zratio1 = zratio1))
})

test_that("bridgeF_nb is increasing in r", {
  expect_lte(bridgeF_nb(r = r[1], zratio1 = zratio1, zratio2 = zratio2), bridgeF_nb(r = r[2], zratio1 = zratio1, zratio2 = zratio2))
})

test_that("bridgeF_nt is increasing in r", {
  expect_lte(bridgeF_nt(r = r[1], zratio1 = zratio1, zratio2 = zratio2), bridgeF_nt(r = r[2], zratio1 = zratio1, zratio2 = zratio2))
})

zratio2 = matrix(c(.2, .8), nrow = 1)

test_that("bridgeF_nn is increasing in r", {
  expect_lte(bridgeF_nn(r = r[1], zratio1 = zratio1, zratio2 = zratio2), bridgeF_nn(r = r[2], zratio1 = zratio1, zratio2 = zratio2))
})
