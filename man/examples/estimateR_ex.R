### Data setting
n <- 1000; p1 <- 15; p2 <- 10 # sample size and dimensions for two datasets.

### Correlation structure within each data set
set.seed(0)
perm1 <- sample(1:(p1 + p2), size = p1 + p2);
Sigma <- autocor(p1 + p2, 0.7)[perm1, perm1]
mu <- rbinom(p1 + p2, 1, 0.5)
# Data generation
simdata <- GenData(n=n, type1 = "ternary", type2 = "continuous", p1 = p1, p2 = p2, copula1 = "exp",
                   copula2 = "cube",  muZ = mu, Sigma = Sigma,
                   c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1), c2 =  NULL)
X1 <- simdata$X1; X2 <- simdata$X2
# Estimate latent correlation matrix with original method
R_nc_org <- estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "continuous",
                 method = "original")$R
# Estimate latent correlation matrix with aprroximation method
R_nc_approx <- estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "continuous",
                    method = "approx")$R

cp1 <- "exp"; cp2 <- "cube"
for (tp1 in c("continuous", "binary", "ternary", "trunc")) {
  for (tp2 in c("continuous", "binary", "ternary", "trunc")) {
    for (md in c("original", "ml", "approx")) {
      c1 <- c2 <- NULL
      if (tp1 == "binary" | tp1 == "trunc") {
        c1 <- rep(1, p1)
      } else if (tp1 == "ternary") {
        c1 <- matrix(rep(1:2, p1), nrow = 2, ncol = p1)
      }
      if (tp2 == "binary" | tp2 == "trunc") {
        c2 <- rep(0, p2)
      } else if (tp2 == "ternary") {
        c2 <- matrix(rep(0:1, p2), nrow = 2, ncol = p2)
      }
      simdata <- GenData(n=n, copula1 = cp1, copula2 = cp2, type1 = tp1, type2 = tp2,
                         muZ = mu, Sigma = Sigma, p1= p1, p2 = p2, c1 = c1, c2 = c2)
      X1 <- simdata$X1; X2 <- simdata$X2
      assign(paste("R", cp1, cp2, tp1, tp2, md, sep = "_"),
             estR(X1 = X1, type1 = tp1, X2 = X2, type2 = tp2, method = md)$R)
      PlotPair(datapair = cbind(c(Sigma), c(get(paste("R", cp1, cp2, tp1, tp2, md, sep = "_")))),
               namepair = c("Sigma", paste("R", cp1, cp2, tp1, tp2, md, sep = "_")),
               title = "Latent correlation (True vs. Estimated)")
    }
  }
}

# Plot ternary/continuous case estimation via original method.
PlotPair(datapair = cbind(c(Sigma), c(R_nc_org)), namepair = c("Sigma", "R_nc_org"),
         title = "Latent correlation (True vs. Estimated)")
# Plot ternary/continuous case estimation via approximation method.
PlotPair(datapair = cbind(c(Sigma), c(R_nc_approx)), namepair = c("Sigma", "R_nc_approx"),
         title = "Latent correlation (True vs. Estimated)")
