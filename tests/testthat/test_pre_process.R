library(latentcor)

test_that("Test non-numeric conversion and types generation", {
  raw_data = data.frame(matrix(NA, 3, 3))
  raw_data$X1 = 1:3
  raw_data$X2 = c("medium", "small", "large")
  raw_data$X3 = 7:9
  pre_process = pre_process(raw_data = raw_data, ordinals = list(NA, c("small", "medium", "large"), NA))
  expect_equal(pre_process$X, data.matrix(data.frame(matrix(c(1, 2, 3, 1, 0, 2, 7, 8, 9), 3, 3))))
  expect_equal(pre_process$types, c("ter", "ter", "ter"))
})
