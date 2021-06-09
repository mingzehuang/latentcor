
library(latentcor)

### Data setting
n <- 1000; p1 <- 15; p2 <- 10 # sample size and dimensions for two datasets.

### Correlation structure within each data set
set.seed(0)
perm1 <- sample(1:(p1 + p2), size = p1 + p2);
Sigma <- autocor(p1 + p2, 0.7)[perm1, perm1]
mu <- rbinom(p1 + p2, 1, 0.5)

XY = mvrnorm(n, mu = mu, Sigma = Sigma)
X = XY[ , 1:p1]; Y = XY[ , (p1 + 1):(p1 + p2)]

test_that("Kendall's matrix is symmetric matrix", {
  expect_equal(Kendall_matrix(X), t(Kendall_matrix(X)))
  expect_equal(Kendall_matrix(X, Y), t(Kendall_matrix(Y, X)))
})

test_that("Kendall's tau is between -1 and 1", {
  expect_true(all(Kendall_matrix(X) <= matrix(1, p1, p1)))
  expect_true(all(Kendall_matrix(X) >= matrix(-1, p1, p1)))
  expect_true(all(c(Kendall_matrix(X, Y)) <= rep(1, p1 * p2)))
  expect_true(all(c(Kendall_matrix(X, Y)) >= rep(-1, p1 * p2)))
})

