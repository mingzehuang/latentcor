# Data generation
X = GenData()$X
# Estimate latent correlation matrix with original method
R_nc_org = latentcor(X = X, types = c("ter", "con"), method = "original")$R
# Estimate latent correlation matrix with aprroximation method
R_nc_approx = latentcor(X = X, types = c("ter", "con"), method = "approx")$R
# Heatmap for latent correlation matrix.
Heatmap_R_nc_approx = latentcor(X = X, types = c("ter", "con"), method = "approx",
                                showplot = TRUE)$plotR
