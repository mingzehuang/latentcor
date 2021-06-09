
library(latentcor)

### Data setting
n <- 1000; p1 <- 15; p2 <- 10 # sample size and dimensions for two datasets.
maxcancor <- 0.9 # true canonical correlation

### Correlation structure within each data set
set.seed(0)
perm1 <- sample(1:(p1 + p2), size = p1 + p2);
Sigma <- autocor(p1 + p2, 0.7)[perm1, perm1]
mu <- rbinom(p1 + p2, 1, 0.5)
# Data generation
simdata <- GenData(n=n, type1 = "continuous", type2 = "binary", p1 = p1, p2 = p2, copula1 = "exp",
                   copula2 = "cube",  muZ = mu, Sigma = Sigma,
                   c1 = NULL, c2 =  rep(0, p2))
X1 <- simdata$X1; X2 <- simdata$X2

test_that("zratio is NULL if type is continuous", {
  expect_true(is.null(zratio(X1, "continuous")))
})

test_that("zratio is the proportion of zeros if type is binary", {
  expect_equal(zratio(X2, "binary"), as.matrix(colMeans(X2 == 0)))
})

# Data generation
simdata <- GenData(n=n, type1 = "trunc", type2 = "ternary", p1 = p1, p2 = p2, copula1 = "exp",
                   copula2 = "cube",  muZ = mu, Sigma = Sigma,
                   c1 = rep(1, p1), c2 =  matrix(rep(0:1, p2), nrow = 2, ncol = p2))
X1 <- simdata$X1; X2 <- simdata$X2

test_that("zratio is the proportion of zeros if type is truncated", {
  expect_equal(zratio(X1, "trunc"), as.matrix(colMeans(X1 == 0)))
})

test_that("zratio is the proportion of zeros and ones if type is ternary", {
  expect_equal(zratio(X2, "ternary"), cbind(colMeans(X2 == 0), 1 - colMeans(X2 == 2)))
})
