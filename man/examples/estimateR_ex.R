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
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "trunc", type2 = "trunc",
                        c1 = rep(1, p1), c2 =  rep(0, p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12 <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_org <- estimateR(X1, type = "trunc", method = "original")$R
R2_org <- estimateR(X2, type = "trunc", method = "original")$R
R12_org <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "trunc", method = "original")$R12
# Plots
df1 <- data.frame(c(Sigma1), c(R1_org))
colnames(df1) = c("Sigma1", "R1_org")
ggplot(df1, aes(Sigma1, R1_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df2 <- data.frame(c(Sigma2), c(R2_org))
colnames(df2) = c("Sigma2", "R2_org")
ggplot(df2, aes(Sigma2, R2_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df12 <- data.frame(c(Sigma12), c(R12_org))
colnames(df12) = c("Sigma12", "R12_org")
ggplot(df12, aes(Sigma12, R12_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
# Estimate latent correlation matrix with faster approximation method
R1_approx <- estimateR(X1, type = "trunc", method = "approx")$R
R2_approx <- estimateR(X2, type = "trunc", method = "approx")$R
R12_approx <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "trunc", method = "approx")$R12
# Plots
df1 <- data.frame(c(Sigma1), c(R1_approx))
colnames(df1) = c("Sigma1", "R1_approx")
ggplot(df1, aes(Sigma1, R1_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df2 <- data.frame(c(Sigma2), c(R2_approx))
colnames(df2) = c("Sigma2", "R2_approx")
ggplot(df2, aes(Sigma2, R2_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df12 <- data.frame(c(Sigma12), c(R12_approx))
colnames(df12) = c("Sigma12", "R12_approx")
ggplot(df12, aes(Sigma12, R12_approx), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "trunc", type2 = "continuous",
                        c1 = rep(1, p1), c2 =  rep(0, p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12 <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_org <- estimateR(X1, type = "trunc", method = "original")$R
R2_org <- estimateR(X2, type = "continuous", method = "original")$R
R12_org <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "continuous", method = "original")$R12
# Plots
df1 <- data.frame(c(Sigma1), c(R1_org))
colnames(df1)=c("Sigma1", "R1_org")
ggplot(df1, aes(Sigma1, R1_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df2 <- data.frame(c(Sigma2), c(R2_org))
colnames(df2)=c("Sigma2", "R2_org")
ggplot(df2, aes(Sigma2, R2_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df12 <- data.frame(c(Sigma12), c(R12_org))
colnames(df12) = c("Sigma12", "R12_org")
ggplot(df12, aes(Sigma12, R12_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
# Estimate latent correlation matrix with faster approximation method
R1_approx <- estimateR(X1, type = "trunc", method = "approx")$R
R2_approx <- estimateR(X2, type = "continuous", method = "approx")$R
R12_approx <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "continuous", method = "approx")$R12

# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "binary", type2 = "binary",
                        c1 = rep(1, p1), c2 =  rep(0, p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12 <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_org <- estimateR(X1, type = "binary", method = "original")$R
R2_org <- estimateR(X2, type = "binary", method = "original")$R
R12_org <- estimateR_mixed(X1, X2, type1 = "binary", type2 = "binary", method = "original")$R12
# Plots
df1 <- data.frame(c(Sigma1), c(R1_org))
colnames(df1)=c("Sigma1", "R1_org")
ggplot(df1, aes(Sigma1, R1_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df2 <- data.frame(c(Sigma2), c(R2_org))
colnames(df2)=c("Sigma2", "R2_org")
ggplot(df2, aes(Sigma2, R2_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df12 <- data.frame(c(Sigma12), c(R12_org))
colnames(df12) = c("Sigma12", "R12_org")
ggplot(df12, aes(Sigma12, R12_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
# Estimate latent correlation matrix with faster approximation method
R1_approx <- estimateR(X1, type = "binary", method = "approx")$R
R2_approx <- estimateR(X2, type = "binary", method = "approx")$R
R12_approx <- estimateR_mixed(X1, X2, type1 = "binary", type2 = "binary", method = "approx")$R12

# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "binary", type2 = "continuous",
                        c1 = rep(1, p1), c2 =  rep(0, p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12 <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_org <- estimateR(X1, type = "binary", method = "original")$R
R2_org <- estimateR(X2, type = "continuous", method = "original")$R
R12_org <- estimateR_mixed(X1, X2, type1 = "binary", type2 = "continuous", method = "original")$R12
# Plots
df1 <- data.frame(c(Sigma1), c(R1_org))
colnames(df1)=c("Sigma1", "R1_org")
ggplot(df1, aes(Sigma1, R1_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df2 <- data.frame(c(Sigma2), c(R2_org))
colnames(df2)=c("Sigma2", "R2_org")
ggplot(df2, aes(Sigma2, R2_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df12 <- data.frame(c(Sigma12), c(R12_org))
colnames(df12) = c("Sigma12", "R12_org")
ggplot(df12, aes(Sigma12, R12_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
# Estimate latent correlation matrix with faster approximation method
R1_approx <- estimateR(X1, type = "binary", method = "approx")$R
R2_approx <- estimateR(X2, type = "continuous", method = "approx")$R
R12_approx <- estimateR_mixed(X1, X2, type1 = "binary", type2 = "continuous", method = "approx")$R12

# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "trunc", type2 = "binary",
                        c1 = rep(1, p1), c2 =  rep(0, p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12 <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_org <- estimateR(X1, type = "trunc", method = "original")$R
R2_org <- estimateR(X2, type = "binary", method = "original")$R
R12_org <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "binary", method = "original")$R12
# Plots
df1 <- data.frame(c(Sigma1), c(R1_org))
colnames(df1)=c("Sigma1", "R1_org")
ggplot(df1, aes(Sigma1, R1_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df2 <- data.frame(c(Sigma2), c(R2_org))
colnames(df2)=c("Sigma2", "R2_org")
ggplot(df2, aes(Sigma2, R2_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df12 <- data.frame(c(Sigma12), c(R12_org))
colnames(df12) = c("Sigma12", "R12_org")
ggplot(df12, aes(Sigma12, R12_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
# Estimate latent correlation matrix with faster approximation method
R1_approx <- estimateR(X1, type = "trunc", method = "approx")$R
R2_approx <- estimateR(X2, type = "binary", method = "approx")$R
R12_approx <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "binary", method = "approx")$R12

# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "ternary", type2 = "ternary",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1), c2 =  matrix(rep(0:1, p2), nrow = 2, ncol = p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12 <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_org <- estimateR(X1, type = "ternary", method = "original")$R
R2_org <- estimateR(X2, type = "ternary", method = "original")$R
R12_org <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "ternary", method = "original")$R12
# Plots
df1 <- data.frame(c(Sigma1), c(R1_org))
colnames(df1)=c("Sigma1", "R1_org")
ggplot(df1, aes(Sigma1, R1_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df2 <- data.frame(c(Sigma2), c(R2_org))
colnames(df2)=c("Sigma2", "R2_org")
ggplot(df2, aes(Sigma2, R2_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df12 <- data.frame(c(Sigma12), c(R12_org))
colnames(df12) = c("Sigma12", "R12_org")
ggplot(df12, aes(Sigma12, R12_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
# Estimate latent correlation matrix with faster approximation method
# R1_approx <- estimateR(X1, type = "ternary", method = "approx")$R
# R2_approx <- estimateR(X2, type = "ternary", method = "approx")$R
# R12_approx <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "ternary", method = "approx")$R12
#
# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "ternary", type2 = "continuous",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1), c2 =  matrix(rep(0:1, p2), nrow = 2, ncol = p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12 <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_org <- estimateR(X1, type = "ternary", method = "original")$R
R2_org <- estimateR(X2, type = "continuous", method = "original")$R
R12_org <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "continuous", method = "original")$R12
# Plots
df1 <- data.frame(c(Sigma1), c(R1_org))
colnames(df1)=c("Sigma1", "R1_org")
ggplot(df1, aes(Sigma1, R1_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df2 <- data.frame(c(Sigma2), c(R2_org))
colnames(df2)=c("Sigma2", "R2_org")
ggplot(df2, aes(Sigma2, R2_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df12 <- data.frame(c(Sigma12), c(R12_org))
colnames(df12) = c("Sigma12", "R12_org")
ggplot(df12, aes(Sigma12, R12_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
# Estimate latent correlation matrix with faster approximation method
# # R1_approx <- estimateR(X1, type = "ternary", method = "approx")$R
R2_approx <- estimateR(X2, type = "continuous", method = "approx")$R
# # R12_approx <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "continuous", method = "approx")$R12

# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "continuous", type2 = "ternary",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1), c2 =  matrix(rep(0:1, p2), nrow = 2, ncol = p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12 <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_org <- estimateR(X1, type = "continuous", method = "original")$R
R2_org <- estimateR(X2, type = "ternary", method = "original")$R
R12_org <- estimateR_mixed(X1, X2, type1 = "continuous", type2 = "ternary", method = "original")$R12
# Plots
df1 <- data.frame(c(Sigma1), c(R1_org))
colnames(df1)=c("Sigma1", "R1_org")
ggplot(df1, aes(Sigma1, R1_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df2 <- data.frame(c(Sigma2), c(R2_org))
colnames(df2)=c("Sigma2", "R2_org")
ggplot(df2, aes(Sigma2, R2_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df12 <- data.frame(c(Sigma12), c(R12_org))
colnames(df12) = c("Sigma12", "R12_org")
ggplot(df12, aes(Sigma12, R12_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
# Estimate latent correlation matrix with faster approximation method
R1_approx <- estimateR(X1, type = "continuous", method = "approx")$R
# # R2_approx <- estimateR(X2, type = "ternary", method = "approx")$R
# # R12_approx <- estimateR_mixed(X1, X2, type1 = "continuous", type2 = "ternary", method = "approx")$R12

# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "ternary", type2 = "binary",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1), c2 = rep(0, p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12 <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_org <- estimateR(X1, type = "ternary", method = "original")$R
R2_org <- estimateR(X2, type = "binary", method = "original")$R
R12_org <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "binary", method = "original")$R12
# Plots
df1 <- data.frame(c(Sigma1), c(R1_org))
colnames(df1)=c("Sigma1", "R1_org")
ggplot(df1, aes(Sigma1, R1_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df2 <- data.frame(c(Sigma2), c(R2_org))
colnames(df2)=c("Sigma2", "R2_org")
ggplot(df2, aes(Sigma2, R2_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df12 <- data.frame(c(Sigma12), c(R12_org))
colnames(df12) = c("Sigma12", "R12_org")
ggplot(df12, aes(Sigma12, R12_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
# Estimate latent correlation matrix with faster approximation method
# # R1_approx <- estimateR(X1, type = "ternary", method = "approx")$R
R2_approx <- estimateR(X2, type = "binary", method = "approx")$R
# # R12_approx <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "binary", method = "approx")$R12

# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "binary", type2 = "ternary",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1), c2 = rep(0, p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12 <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_org <- estimateR(X1, type = "binary", method = "original")$R
R2_org <- estimateR(X2, type = "ternary", method = "original")$R
R12_org <- estimateR_mixed(X1, X2, type1 = "binary", type2 = "ternary", method = "original")$R12
# Plots
df1 <- data.frame(c(Sigma1), c(R1_org))
colnames(df1)=c("Sigma1", "R1_org")
ggplot(df1, aes(Sigma1, R1_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df2 <- data.frame(c(Sigma2), c(R2_org))
colnames(df2)=c("Sigma2", "R2_org")
ggplot(df2, aes(Sigma2, R2_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df12 <- data.frame(c(Sigma12), c(R12_org))
colnames(df12) = c("Sigma12", "R12_org")
ggplot(df12, aes(Sigma12, R12_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
# Estimate latent correlation matrix with faster approximation method
R1_approx <- estimateR(X1, type = "binary", method = "approx")$R
# # R2_approx <- estimateR(X2, type = "ternary", method = "approx")$R
# # R12_approx <- estimateR_mixed(X1, X2, type1 = "binary", type2 = "ternary", method = "approx")$R12

# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "ternary", type2 = "trunc",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1), c2 = rep(0, p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12 <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_org <- estimateR(X1, type = "ternary", method = "original")$R
R2_org <- estimateR(X2, type = "trunc", method = "original")$R
# # R12_org <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "trunc", method = "original")$R12
# Plots
df1 <- data.frame(c(Sigma1), c(R1_org))
colnames(df1)=c("Sigma1", "R1_org")
ggplot(df1, aes(Sigma1, R1_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df2 <- data.frame(c(Sigma2), c(R2_org))
colnames(df2)=c("Sigma2", "R2_org")
ggplot(df2, aes(Sigma2, R2_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df12 <- data.frame(c(Sigma12), c(R12_org))
colnames(df12) = c("Sigma12", "R12_org")
ggplot(df12, aes(Sigma12, R12_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
# Estimate latent correlation matrix with faster approximation method
# # R1_approx <- estimateR(X1, type = "ternary", method = "approx")$R
R2_approx <- estimateR(X2, type = "trunc", method = "approx")$R
# # R12_approx <- estimateR_mixed(X1, X2, type1 = "ternary", type2 = "trunc", method = "approx")$R12

# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "trunc", type2 = "ternary",
                        c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1), c2 = rep(0, p2)
)
X1 <- simdata$X1
X2 <- simdata$X2
Sigma12 <- simdata$Sigma12
# Estimate latent correlation matrix with original method
R1_org <- estimateR(X1, type = "trunc", method = "original")$R
R2_org <- estimateR(X2, type = "ternary", method = "original")$R
# R12_org <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "ternary", method = "original")$R12
# Plots
df1 <- data.frame(c(Sigma1), c(R1_org))
colnames(df1)=c("Sigma1", "R1_org")
ggplot(df1, aes(Sigma1, R1_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df2 <- data.frame(c(Sigma2), c(R2_org))
colnames(df2)=c("Sigma2", "R2_org")
ggplot(df2, aes(Sigma2, R2_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
df12 <- data.frame(c(Sigma12), c(R12_org))
colnames(df12) = c("Sigma12", "R12_org")
ggplot(df12, aes(Sigma12, R12_org), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
# # Estimate latent correlation matrix with faster approximation method
R1_approx <- estimateR(X1, type = "trunc", method = "approx")$R
# R2_approx <- estimateR(X2, type = "ternary", method = "approx")$R
# R12_approx <- estimateR_mixed(X1, X2, type1 = "trunc", type2 = "ternary", method = "approx")$R12

### Check the range of truncation levels of variables
range(colMeans(X1 == 0))
range(colMeans(X2 == 0))



