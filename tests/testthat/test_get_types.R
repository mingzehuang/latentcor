library(latentcor)

test_that("Test non-numeric conversion and types generation", {
  X = gen_data(types = c("ter", "con"))$X
  expect_equal(get_types(X = X), c("ter", "con"))
})
