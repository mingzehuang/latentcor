
library(latentcor)

r = runif(1, -1, 1); zratio1 = .5; zratio2 = .5

test_that("Bridge selection should be consistent with each bridge function.", {
  expect_equal(bridge(type1 = "binary", type2 = "continuous", r = r, zratio1 = zratio1), bridgeF_bc(r = r, zratio1 = zratio1))
  expect_equal(bridge(type1 = "binary", type2 = "binary", r = r, zratio1 = zratio1, zratio2 = zratio2), bridgeF_bb(r = r, zratio1 = zratio1, zratio2 = zratio2))
  expect_equal(bridge(type1 = "trunc", type2 = "continuous", r = r, zratio1 = zratio1), bridgeF_tc(r = r, zratio1 = zratio1))
  expect_equal(bridge(type1 = "trunc", type2 = "binary", r = r, zratio1 = zratio1, zratio2 = zratio2), bridgeF_tb(r = r, zratio1 = zratio1, zratio2 = zratio2))
  expect_equal(bridge(type1 = "trunc", type2 = "trunc", r = r, zratio1 = zratio1, zratio2 = zratio2), bridgeF_tt(r = r, zratio1 = zratio1, zratio2 = zratio2))
})

zratio1 = c(.3, .7)
test_that("Bridge selection should be consistent with each bridge function.", {
  expect_equal(bridge(type1 = "ternary", type2 = "continuous", r = r, zratio1 = zratio1), bridgeF_nc(r = r, zratio1 = zratio1))
  expect_equal(bridge(type1 = "ternary", type2 = "binary", r = r, zratio1 = zratio1, zratio2 = zratio2), bridgeF_nb(r = r, zratio1 = zratio1, zratio2 = zratio2))
  expect_equal(bridge(type1 = "ternary", type2 = "trunc", r = r, zratio1 = zratio1, zratio2 = zratio2), bridgeF_nt(r = r, zratio1 = zratio1, zratio2 = zratio2))
})

zratio2 = c(.2, .8)
test_that("Bridge selection should be consistent with each bridge function.", {
  expect_equal(bridge(type1 = "ternary", type2 = "ternary", r = r, zratio1 = zratio1, zratio2 = zratio2), bridgeF_nn(r = r, zratio1 = zratio1, zratio2 = zratio2))
})
