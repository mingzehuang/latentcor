context("Test KendallTau")
library(latentcor)

test_that("Kendall' tau is symmetric for two inputs", {
  X1 = rnorm(10); X2 = rnorm(10)
  expect_equal(KendallTau(X1, X2), KendallTau(X2, X1))
})
#> Test passed 😸

#> Test passed 🌈

#> Test passed 🎉
