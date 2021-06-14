### Data setting
n = 1000; p1 = 15; p2 = 10 # sample size and dimensions for two datasets.
rho = .9
# Data generation
simdata = GenData(n=n, type1 = "ternary", type2 = "continuous", p1 = p1, p2 = p2, rho = rho,
                  copula1 = "cube", copula2 = "cube", c1 = c(0, 1), c2 =  NULL)
Sigma = simdata$Sigma; X1 = simdata$X1; X2 = simdata$X2
# Estimate latent correlation matrix with original method
R_nc_org = estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "continuous", method = "original")$R
# Estimate latent correlation matrix with aprroximation method
R_nc_approx = estR(X1 = X1, type1 = "ternary", X2 = X2, type2 = "continuous", method = "approx")$R

