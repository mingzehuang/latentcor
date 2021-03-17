
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
                        type1 = "trunc", type2 = "trunc",
                        c1 = rep(1, p1), c2 =  rep(0, p2)
)
X1 <- simdata$X1; X2 <- simdata$X2; Sigma12_tt <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_tt_org <- estimateR(X1, type = "trunc", method = "original")$R
R2_tt_org <- estimateR(X2, type = "trunc", method = "original")$R
R12_tt_org <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "trunc",
                              method = "original")$R12;
for (cp1 in c("exp", "cube")) {
  for (cp2 in c("exp", "cube")) {
    for (tp1 in c("continuous", "binary", "trunc")) {
      for (tp2 in c("continuous", "binary", "trunc")) {
        for (md in c("original", "approx")) {
          simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                                    maxcancor = maxcancor,
                                    Sigma1 = Sigma1, Sigma2 = Sigma2,
                                    copula1 = cp1, copula2 = cp2,
                                    muZ = mu,
                                    type1 = tp1, type2 = tp2,
                                    c1 = rep(1, p1), c2 =  rep(0, p2))
          assign(paste("R1", cp1, cp2, tp1, tp2, md, sep = "_"),
                 estimateR(X1, type = tp1, method = md)$R)
          assign(paste("R2", cp1, cp2, tp1, tp2, md, sep = "_"),
                 estimateR(X2, type = tp2, method = md)$R)
          assign(paste("R12", cp1, cp2, tp1, tp2, md, sep = "_"),
                 estimateR_mixed(X1, X2, type1 = tp1, type2 = tp2,
                                 method = md)$R12)
        }
      }
    }
  }
}
# Estimate latent correlation matrix with faster approximation method
R1_tt_approx <- estimateR(X1, type = "trunc", method = "approx")$R
R2_tt_approx <- estimateR(X2, type = "trunc", method = "approx")$R
R12_tt_approx <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "trunc",
                                 method = "approx")$R12
# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                        maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "trunc", type2 = "continuous",
                        c1 = rep(1, p1), c2 =  rep(0, p2)
)
X1 <- simdata$X1; X2 <- simdata$X2; Sigma12_tc <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_tc_org <- estimateR(X1, type = "trunc", method = "original")$R
R2_tc_org <- estimateR(X2, type = "continuous", method = "original")$R
R12_tc_org <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "continuous",
                              method = "original")$R12
# Estimate latent correlation matrix with faster approximation method
R1_tc_approx <- estimateR(X1, type = "trunc", method = "approx")$R
R2_tc_approx <- estimateR(X2, type = "continuous", method = "approx")$R
R12_tc_approx <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "continuous",
                                 method = "approx")$R12
# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                        maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "binary", type2 = "binary",
                        c1 = rep(1, p1), c2 =  rep(0, p2)
)
X1 <- simdata$X1; X2 <- simdata$X2; Sigma12_bb <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_bb_org <- estimateR(X1, type = "binary", method = "original")$R
R2_bb_org <- estimateR(X2, type = "binary", method = "original")$R
R12_bb_org <- estimateR_mixed(X1, X2, type1 = "binary", type2 = "binary",
                              method = "original")$R12
# Estimate latent correlation matrix with faster approximation method
R1_bb_approx <- estimateR(X1, type = "binary", method = "approx")$R
R2_bb_approx <- estimateR(X2, type = "binary", method = "approx")$R
R12_bb_approx <- estimateR_mixed(X1, X2, type1 = "binary", type2 = "binary",
                                 method = "approx")$R12
# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                        maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "binary", type2 = "continuous",
                        c1 = rep(1, p1), c2 =  rep(0, p2)
)
X1 <- simdata$X1; X2 <- simdata$X2; Sigma12_bc <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_bc_org <- estimateR(X1, type = "binary", method = "original")$R
R2_bc_org <- estimateR(X2, type = "continuous", method = "original")$R
R12_bc_org <- estimateR_mixed(X1, X2, type1 = "binary", type2 = "continuous",
                              method = "original")$R12
# Estimate latent correlation matrix with faster approximation method
R1_bc_approx <- estimateR(X1, type = "binary", method = "approx")$R
R2_bc_approx <- estimateR(X2, type = "continuous", method = "approx")$R
R12_bc_approx <- estimateR_mixed(X1, X2, type1 = "binary", type2 = "continuous",
                                 method = "approx")$R12
# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                        maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "trunc", type2 = "binary",
                        c1 = rep(1, p1), c2 =  rep(0, p2)
)
X1 <- simdata$X1; X2 <- simdata$X2; Sigma12_tb <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_tb_org <- estimateR(X1, type = "trunc", method = "original")$R
R2_tb_org <- estimateR(X2, type = "binary", method = "original")$R
R12_tb_org <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "binary",
                              method = "original")$R12
# Estimate latent correlation matrix with faster approximation method
R1_tb_approx <- estimateR(X1, type = "trunc", method = "approx")$R
R2_tb_approx <- estimateR(X2, type = "binary", method = "approx")$R
R12_tb_approx <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "binary",
                                 method = "approx")$R12
# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                        maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "ternary", type2 = "ternary",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1),
                        c2 =  matrix(rep(0:1, p2), nrow = 2, ncol = p2)
)
X1 <- simdata$X1; X2 <- simdata$X2; Sigma12_nn <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_nn_org <- estimateR(X1, type = "ternary", method = "original")$R
R2_nn_org <- estimateR(X2, type = "ternary", method = "original")$R
R12_nn_org <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "ternary",
                              method = "original")$R12
# Estimate latent correlation matrix with faster approximation method
R1_nn_approx <- estimateR(X1, type = "ternary", method = "approx")$R
R2_nn_approx <- estimateR(X2, type = "ternary", method = "approx")$R
R12_nn_approx <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "ternary",
                                 method = "approx")$R12
# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                        maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "ternary", type2 = "continuous",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1),
                        c2 = matrix(rep(0:1, p2), nrow = 2, ncol = p2)
)
X1 <- simdata$X1; X2 <- simdata$X2; Sigma12_nc <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_nc_org <- estimateR(X1, type = "ternary", method = "original")$R
R2_nc_org <- estimateR(X2, type = "continuous", method = "original")$R
R12_nc_org <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "continuous",
                              method = "original")$R12
# Estimate latent correlation matrix with faster approximation method
R1_nc_approx <- estimateR(X1, type = "ternary", method = "approx")$R
R2_nc_approx <- estimateR(X2, type = "continuous", method = "approx")$R
R12_nc_approx <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "continuous",
                                 method = "approx")$R12
# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                        maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "continuous", type2 = "ternary",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1),
                        c2 = matrix(rep(0:1, p2), nrow = 2, ncol = p2)
)
X1 <- simdata$X1; X2 <- simdata$X2; Sigma12_cn <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_cn_org <- estimateR(X1, type = "continuous", method = "original")$R
R2_cn_org <- estimateR(X2, type = "ternary", method = "original")$R
R12_cn_org <- estimateR_mixed(X1, X2, type1 = "continuous", type2 = "ternary",
                              method = "original")$R12
# Estimate latent correlation matrix with faster approximation method
R1_cn_approx <- estimateR(X1, type = "continuous", method = "approx")$R
R2_cn_approx <- estimateR(X2, type = "ternary", method = "approx")$R
R12_cn_approx <- estimateR_mixed(X1, X2, type1 = "continuous", type2 = "ternary",
                                 method = "approx")$R12
# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                        maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "ternary", type2 = "binary",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1),
                        c2 = rep(0, p2)
)
X1 <- simdata$X1; X2 <- simdata$X2; Sigma12_nb <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_nb_org <- estimateR(X1, type = "ternary", method = "original")$R
R2_nb_org <- estimateR(X2, type = "binary", method = "original")$R
R12_nb_org <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "binary",
                              method = "original")$R12
# Estimate latent correlation matrix with faster approximation method
R1_nb_approx <- estimateR(X1, type = "ternary", method = "approx")$R
R2_nb_approx <- estimateR(X2, type = "binary", method = "approx")$R
R12_nb_approx <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "binary",
                                 method = "approx")$R12
# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                        maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "binary", type2 = "ternary",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1),
                        c2 = rep(0, p2)
)
X1 <- simdata$X1; X2 <- simdata$X2; Sigma12_bn <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_bn_org <- estimateR(X1, type = "binary", method = "original")$R
R2_bn_org <- estimateR(X2, type = "ternary", method = "original")$R
R12_bn_org <- estimateR_mixed(X1, X2, type1 = "binary", type2 = "ternary",
                              method = "original")$R12
# Estimate latent correlation matrix with faster approximation method
R1_bn_approx <- estimateR(X1, type = "binary", method = "approx")$R
# R2_bn_approx <- estimateR(X2, type = "ternary", method = "approx")$R
R12_bn_approx <- estimateR_mixed(X1, X2, type1 = "binary", type2 = "ternary",
                                 method = "approx")$R12
# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                        maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "ternary", type2 = "trunc",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1),
                        c2 = rep(0, p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12_nt <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_nt_org <- estimateR(X1, type = "ternary", method = "original")$R
R2_nt_org <- estimateR(X2, type = "trunc", method = "original")$R
# R12_nt_org <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "trunc",
#                               method = "original")$R12
# Estimate latent correlation matrix with faster approximation method
# R1_nt_approx <- estimateR(X1, type = "ternary", method = "approx")$R
R2_nt_approx <- estimateR(X2, type = "trunc", method = "approx")$R
# R12_nt_approx <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "trunc",
#                                  method = "approx")$R12
# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2,
                        maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "trunc", type2 = "ternary",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1),
                        c2 = rep(0, p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12_tn <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_tn_org <- estimateR(X1, type = "trunc", method = "original")$R
R2_tn_org <- estimateR(X2, type = "ternary", method = "original")$R
# R12_tn_org <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "ternary",
#                               method = "original")$R12
# Estimate latent correlation matrix with faster approximation method
R1_tn_approx <- estimateR(X1, type = "trunc", method = "approx")$R
# R2_tn_approx <- estimateR(X2, type = "ternary", method = "approx")$R
# R12_tn_approx <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "ternary",
#                                  method = "approx")$R12
### Check the range of truncation levels of variables
range(colMeans(X1 == 0))
range(colMeans(X2 == 0))

# Plot true latent correlation vs. estimated latent correlation
Plot <- PlotCompare(list(cbind(c(Sigma1), c(R1_tt_org)),
                         cbind(c(Sigma2), c(R2_tt_org)),
                         cbind(c(Sigma12_tt), c(R12_tt_org)),
                         cbind(c(Sigma1), c(R1_tt_approx)),
                         cbind(c(Sigma2), c(R2_tt_approx)),
                         cbind(c(Sigma12_tt), c(R12_tt_approx)),
                         cbind(c(Sigma1), c(R1_tc_org)),
                         cbind(c(Sigma2), c(R2_tc_org)),
                         cbind(c(Sigma12_tc), c(R12_tc_org)),
                         cbind(c(Sigma1), c(R1_tc_approx)),
                         cbind(c(Sigma2), c(R2_tc_approx)),
                         cbind(c(Sigma12_tc), c(R12_tc_approx)),
                         cbind(c(Sigma1), c(R1_bb_org)),
                         cbind(c(Sigma2), c(R2_bb_org)),
                         cbind(c(Sigma12_bb), c(R12_bb_org)),
                         cbind(c(Sigma1), c(R1_bb_approx)),
                         cbind(c(Sigma2), c(R2_bb_approx)),
                         cbind(c(Sigma12_bb), c(R12_bb_approx)),
                         cbind(c(Sigma1), c(R1_bc_org)),
                         cbind(c(Sigma2), c(R2_bc_org)),
                         cbind(c(Sigma12_bc), c(R12_bc_org)),
                         cbind(c(Sigma1), c(R1_bc_approx)),
                         cbind(c(Sigma2), c(R2_bc_approx)),
                         cbind(c(Sigma12_bc), c(R12_bc_approx)),
                         cbind(c(Sigma1), c(R1_tb_org)),
                         cbind(c(Sigma2), c(R2_tb_org)),
                         cbind(c(Sigma12_tb), c(R12_tb_org)),
                         cbind(c(Sigma1), c(R1_tb_approx)),
                         cbind(c(Sigma2), c(R2_tb_approx)),
                         cbind(c(Sigma12_tb), c(R12_tb_approx)),
                         cbind(c(Sigma1), c(R1_nn_org)),
                         cbind(c(Sigma2), c(R2_nn_org)),
                         cbind(c(Sigma12_nn), c(R12_nn_org)),
                         cbind(c(Sigma1), c(R1_nn_approx)),
                         cbind(c(Sigma2), c(R2_nn_approx)),
                         cbind(c(Sigma12_nn), c(R12_nn_approx)),
                         cbind(c(Sigma1), c(R1_nc_org)),
                         cbind(c(Sigma2), c(R2_nc_org)),
                         cbind(c(Sigma12_nc), c(R12_nc_org)),
                         cbind(c(Sigma1), c(R1_nc_approx)),
                         cbind(c(Sigma2), c(R2_nc_approx)),
                         cbind(c(Sigma12_nc), c(R12_nc_approx)),
                         cbind(c(Sigma1), c(R1_cn_org)),
                         cbind(c(Sigma2), c(R2_cn_org)),
                         cbind(c(Sigma12_cn), c(R12_cn_org)),
                         cbind(c(Sigma1), c(R1_cn_approx)),
                         cbind(c(Sigma2), c(R2_cn_approx)),
                         cbind(c(Sigma12_cn), c(R12_cn_approx)),
                         cbind(c(Sigma1), c(R1_nb_org)),
                         cbind(c(Sigma2), c(R2_nb_org)),
                         cbind(c(Sigma12_nb), c(R12_nb_org)),
                         cbind(c(Sigma1), c(R1_nb_approx)),
                         cbind(c(Sigma2), c(R2_nb_approx)),
                         cbind(c(Sigma12_nb), c(R12_nb_approx)),
                         cbind(c(Sigma1), c(R1_bn_org)),
                         cbind(c(Sigma2), c(R2_bn_org)),
                         cbind(c(Sigma12_bn), c(R12_bn_org)),
                         cbind(c(Sigma1), c(R1_bn_approx)),
                         cbind(c(Sigma12_bn), c(R12_bn_approx)),
                         cbind(c(Sigma1), c(R1_nt_org)),
                         cbind(c(Sigma2), c(R2_nt_org)),
                         cbind(c(Sigma2), c(R2_nt_approx)),
                         cbind(c(Sigma1), c(R1_tn_org)),
                         cbind(c(Sigma2), c(R2_tn_org)),
                         cbind(c(Sigma1), c(R1_tn_approx))
),
list(c("Sigma1", "R1_tt_org"),
     c("Sigma2", "R2_tt_org"),
     c("Sigma12_tt", "R12_tt_org"),
     c("Sigma1", "R1_tt_approx"),
     c("Sigma2", "R2_tt_approx"),
     c("Sigma12_tt", "R12_tt_approx"),
     c("Sigma1", "R1_tc_org"),
     c("Sigma2", "R2_tc_org"),
     c("Sigma12_tc", "R12_tc_org"),
     c("Sigma1", "R1_tc_approx"),
     c("Sigma2", "R2_tc_approx"),
     c("Sigma12_tc", "R12_tc_approx"),
     c("Sigma1", "R1_bb_org"),
     c("Sigma2", "R2_bb_org"),
     c("Sigma12_bb", "R12_bb_org"),
     c("Sigma1", "R1_bb_approx"),
     c("Sigma2", "R2_bb_approx"),
     c("Sigma12_bb", "R12_bb_approx"),
     c("Sigma1", "R1_bc_org"),
     c("Sigma2", "R2_bc_org"),
     c("Sigma12_bc", "R12_bc_org"),
     c("Sigma1", "R1_bc_approx"),
     c("Sigma2", "R2_bc_approx"),
     c("Sigma12_bc", "R12_bc_approx"),
     c("Sigma1", "R1_tb_org"),
     c("Sigma2", "R2_tb_org"),
     c("Sigma12_tb", "R12_tb_org"),
     c("Sigma1", "R1_tb_approx"),
     c("Sigma2", "R2_tb_approx"),
     c("Sigma12_tb", "R12_tb_approx"),
     c("Sigma1", "R1_nn_org"),
     c("Sigma2", "R2_nn_org"),
     c("Sigma12_nn", "R12_nn_org"),
     c("Sigma1", "R1_nn_approx"),
     c("Sigma2", "R2_nn_approx"),
     c("Sigma12_nn", "R12_nn_approx"),
     c("Sigma1", "R1_nc_org"),
     c("Sigma2", "R2_nc_org"),
     c("Sigma12_nc", "R12_nc_org"),
     c("Sigma1", "R1_nc_approx"),
     c("Sigma2", "R2_nc_approx"),
     c("Sigma12_nc", "R12_nc_approx"),
     c("Sigma1", "R1_cn_org"),
     c("Sigma2", "R2_cn_org"),
     c("Sigma12_cn", "R12_cn_org"),
     c("Sigma1", "R1_cn_approx"),
     c("Sigma2", "R2_cn_approx"),
     c("Sigma12_cn", "R12_cn_approx"),
     c("Sigma1", "R1_nb_org"),
     c("Sigma2", "R2_nb_org"),
     c("Sigma12_nb", "R12_nb_org"),
     c("Sigma1", "R1_nb_approx"),
     c("Sigma2", "R2_nb_approx"),
     c("Sigma12_nb", "R12_nb_approx"),
     c("Sigma1", "R1_bn_org"),
     c("Sigma2", "R2_bn_org"),
     c("Sigma12_bn", "R12_bn_org"),
     c("Sigma1", "R1_bn_approx"),
     c("Sigma12_bn", "R12_bn_approx"),
     c("Sigma1", "R1_nt_org"),
     c("Sigma2", "R2_nt_org"),
     c("Sigma2", "R2_nt_approx"),
     c("Sigma1", "R1_tn_org"),
     c("Sigma2", "R2_tn_org"),
     c("Sigma1", "R1_tn_approx")
)
)
