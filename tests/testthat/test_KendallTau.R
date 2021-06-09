
library(latentcor)

X1 = rnorm(10); X2 = rnorm(10)

test_that("Kendall's tau is symmetric for two inputs", {
  expect_equal(KendallTau(X1, X2), KendallTau(X2, X1))
})

test_that("Kendall's tau is between -1 and 1", {
  expect_lte(KendallTau(X1, X2), 1)
  expect_gte(KendallTau(X1, X2), -1)
})

test_that("Kendall's tau captures the order of pairs", {
  X1_asd = sort(X1); X2_asd = sort(X2)
  expect_lte(KendallTau(X1, X2), KendallTau(X1_asd, X2_asd))
})
