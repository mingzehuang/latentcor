
library(ggplot2)
### Data setting
n <- 1000; p1 <- 15; p2 <- 10 # sample size and dimensions for two datasets.
maxcancor <- 0.9 # true canonical correlation

### Correlation structure within each data set
set.seed(0)
perm1 <- sample(1:p1, size = p1);
Sigma1 <- autocor(p1, 0.7)[perm1, perm1]
blockind <- sample(1:3, size = p2, replace = TRUE);
Sigma2 <- blockcor(blockind, 0.7)
mu <- rbinom(p1+p2, 1, 0.5)

### true variable indices for each dataset
trueidx1 <- c(rep(1, 3), rep(0, p1-3))
trueidx2 <- c(rep(1, 2), rep(0, p2-2))

# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                        maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "binary", type2 = "continuous",
                        c1 = rep(1, p1), c2 =  NULL
)
X1 <- simdata$X1; X2 <- simdata$X2; Sigma12_tt <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_tt_org <- estimateR(X1, type = "binary", method = "original")$R
R2_tt_org <- estimateR(X2, type = "continuous", method = "original")$R
R12_tt_org <- estimateR_mixed(X1, X2, type1 = "binary", type2 = "continuous",
                              method = "original")$R12;
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
      simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                                    maxcancor = maxcancor,
                                    Sigma1 = Sigma1, Sigma2 = Sigma2,
                                    copula1 = cp1, copula2 = cp2,
                                    muZ = mu,
                                    type1 = tp1, type2 = tp2,
                                    c1 = c1, c2 = c2)
      X1 <- simdata$X1; X2 <- simdata$X2; Sigma12 <- simdata$Sigma12
      assign(paste("R1", cp1, cp2, tp1, tp2, md, sep = "_"),
             estimateR(X1, type = tp1, method = md)$R)
      assign(paste("R2", cp1, cp2, tp1, tp2, md, sep = "_"),
             estimateR(X2, type = tp2, method = md)$R)
      assign(paste("R12", cp1, cp2, tp1, tp2, md, sep = "_"),
             estimateR_mixed(X1, X2, type1 = tp1, type2 = tp2, method = md)$R12)
      PlotCompare(list(cbind(c(Sigma1), c(get(paste("R1", cp1, cp2, tp1, tp2, md, sep = "_")))),
                       cbind(c(Sigma2), c(get(paste("R2", cp1, cp2, tp1, tp2, md, sep = "_")))),
                       cbind(c(Sigma12), c(get(paste("R12", cp1, cp2, tp1, tp2, md, sep = "_"))))
                        ),
                  list(c("Sigma1", paste("R1", cp1, cp2, tp1, tp2, md, sep = "_")),
                       c("Sigma2", paste("R2", cp1, cp2, tp1, tp2, md, sep = "_")),
                       c("Sigma12", paste("R12", cp1, cp2, tp1, tp2, md, sep = "_"))),
                  "Latent correlation (True vs. Estimated)")
    }
  }
}
#
# #Data generation
# simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
#                         maxcancor = maxcancor,
#                         Sigma1 = Sigma1, Sigma2 = Sigma2,
#                         copula1 = "exp", copula2 = "cube",
#                         muZ = mu,
#                         type1 = "continuous", type2 = "dtrunc",
#                         c1 =  NULL,
#                         c2 <- matrix(rep(0:1, p2), nrow = 2, ncol = p2)
#
# )
# X1 <- simdata$X1; X2 <- simdata$X2; Sigma12_dc <- simdata$Sigma12
# # # Estimate latent correlation matrix with original method
# # R1_dc_org <- estimateR(X1, type = "continuous", method = "original")$R
# # # R2_dc_org <- estimateR(X2, type = "dtrunc", method = "original")$R
# R12_dc_org <- estimateR_mixed(X1, X2, type1 = "continuous", type2 = "dtrunc",
#                               method = "original")$R12
# # Estimate latent correlation matrix with faster approximation method
# R1_dc_approx <- estimateR(X1, type = "dtrunc", method = "approx")$R
# R2_dc_approx <- estimateR(X2, type = "ternary", method = "approx")$R
# R12_dc_approx <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "ternary",
#                                  method = "approx")$R12

# PlotCompare(list(cbind(c(Sigma1), c(R1_dc_org)),
#                  cbind(c(Sigma2), c(R2_dc_org)),
#                  cbind(c(Sigma12), c(R12_dc_org))
# ),
# list(c("Sigma1", "R1_dc_org"),
#      c("Sigma2", "R2_dc_org"),
#      c("Sigma12", "R12_dc_org")),
# "Latent correlation (True vs. Estimated)")

# PlotCompare(list(cbind(c(Sigma1), c(R1_nt_ml)),
#                  cbind(c(Sigma2), c(R2_nt_ml)),
#                  cbind(c(Sigma12), c(R12_nt_ml))
# ),
# list(c("Sigma1", "R1_nt_ml"),
#      c("Sigma2", "R2_nt_ml"),
#      c("Sigma12", "R12_nt_ml")),
# "Latent correlation (True vs. Estimated)")
#
# PlotCompare(list(cbind(c(Sigma1), c(R1_nt_approx)),
#                  cbind(c(Sigma2), c(R2_nt_approx)),
#                  cbind(c(Sigma12), c(R12_nt_approx))
# ),
# list(c("Sigma1", "R1_nt_approx"),
#      c("Sigma2", "R2_nt_approx"),
#      c("Sigma12", "R12_nt_approx")),
# "Latent correlation (True vs. Estimated)")

### Check the range of truncation levels of variables
range(colMeans(X1 == 0))
range(colMeans(X2 == 0))

