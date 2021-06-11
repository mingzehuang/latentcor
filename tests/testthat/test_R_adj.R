
library(latentcor)

R = autocor(3, 0.7)

test_that("Adjust nearPD won't destroy symmetry.", {
  expect_equal(R_adj(R = R, use.nearPD = TRUE, verbose = FALSE, nu = .5),
               t(R_adj(R = R, use.nearPD = TRUE, verbose = FALSE, nu = .5)))
})

R = diag(1, 3); R[3, 3] = 0

test_that("Adjust nearPD won't destroy symmetry.", {
  expect_equal(R_adj(R = R, use.nearPD = TRUE, verbose = FALSE, nu = .5),
               t(R_adj(R = R, use.nearPD = TRUE, verbose = FALSE, nu = .5)))
})
