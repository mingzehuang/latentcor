# Data generation
X = GenData()
# Estimate latent correlation matrix with original method
R_nc_org = estR(X = X, types = c("tru", "ter"), method = "original")$R
# Estimate latent correlation matrix with aprroximation method
R_nc_approx = estR(X = X, types = c("tru", "ter"), method = "approx")$R

# Heatmap for latent correlation matrix.
Heatmap_R_nc_approx = estR(X = X, types = c("tru", "ter"), method = "approx", corplot = TRUE)$plotR
