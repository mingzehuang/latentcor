library(latentcor)
for (comb in latentcor:::combs) {
  assign(paste("ipol", comb, sep = "_"), interpolation(evalfun = latentcor:::evalfun, grid_list = get(paste("latentcor:::grid_list", comb, sep = "_")), comb = comb, tol = 1e-8, ratio = .9)$interpolant)
}
save(list = c(paste("ipol", combs, sep = "_")), file = "interpolation.rda", compress = "xz")
