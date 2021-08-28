library(latentcor)

test_that("Always recode binary/ternary/truncated", {
  types = c("con", "bin", "tru", "ter")
  X = gen_data(types = types)$X; X[sample(length(X), 10)] = NA
  expect_equal(encodeX(X = cbind(X[ , 1], X[ , 2:4] - 5), types = types)$X, X)
})
