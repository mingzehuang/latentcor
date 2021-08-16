## The following script generate interpolant on my HPRC. Users should modify source directory and cores to replicate it.
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/internal.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/interpolation.R")

library(MASS)
library(microbenchmark)
library(foreach)
library(Matrix)
library(chebpol)
library(pcaPP)
library(doRNG)
library(doFuture)

evalfun = function(grid_input, comb, tol, ratio) {
  if (comb == "10" | comb == "20") {
    zratio1 = grid_input[2]; zratio2 = NA
  } else if (comb == "11" | comb == "21" | comb == "22") {
    zratio1 = grid_input[2]; zratio2 = grid_input[3]
  } else if (comb == "30") {
    zratio1 = matrix(c(grid_input[2] * grid_input[3], grid_input[3]), ncol = 1); zratio2 = NA
  } else if (comb == "31" | comb == "32") {
    zratio1 = matrix(c(grid_input[2] * grid_input[3], grid_input[3]), ncol = 1); zratio2 = grid_input[4]
  } else if (comb == "33") {
    zratio1 = matrix(c(grid_input[2] * grid_input[3], grid_input[3]), ncol = 1); zratio2 = matrix(c(grid_input[4:5]), ncol = 1)
  }
  K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
  out = r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = tol, ratio = ratio)
  return(out)
}

grid_list_10 = list(round(pnorm(seq(-1.2, 1.2, by =.06), sd = .5), 6) * 2 - 1, round(pnorm(seq(-1.2, 1.2, by =.06), sd = .5), 6))
grid_list_11 = list(round(pnorm(seq(-1.2, 1.2, by =.06), sd = .5), 6) * 2 - 1, round(pnorm(seq(-1.2, 1.2, by =.06), sd = .5), 6), round(pnorm(seq(-1.2, 1.2, by =.06), sd = .5), 6))
grid_list_20 = list(round(pnorm(seq(-1.2, 1.2, by =.06), sd = .5), 6) * 2 - 1, round(pnorm(seq(.1, 2.5, by = .05)), 6) * 2 - 1)
grid_list_21 = list(round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6) * 2 - 1, round(pnorm(seq(-1.2, 1.2, by =.06), sd = .5), 6), round(pnorm(seq(-1.2, 1.2, by =.06), sd = .5), 6))
grid_list_22 = list(round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6) * 2 - 1, round(pnorm(seq(.1, 2.5, by = .05)), 6) * 2 - 1, round(pnorm(seq(.1, 2.5, by = .05)), 6) * 2 - 1)
grid_list_30 = list(round(pnorm(seq(-1.2, 1.2, by =.06), sd = .5), 6) * 2 - 1, round(pnorm(seq(-2.1, 2.1, by =.1)), 6), round(pnorm(seq(-2.1, 2.1, by =.1)), 6))
grid_list_31 = list(round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6) * 2 - 1, round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6), round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6), round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6))
grid_list_32 = list(round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6) * 2 - 1, round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6), round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6), round(pnorm(seq(.1, 2.5, by = .1)), 6) * 2 - 1)
grid_list_33 = list(round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6) * 2 - 1, round(pnorm(seq(-1.8, 1.8, by =.3), sd = .8), 6), round(pnorm(seq(-1.8, 1.8, by =.3), sd = .8), 6), round(pnorm(seq(-1.8, 1.8, by =.3), sd = .8), 6), round(pnorm(seq(-1.8, 1.8, by =.3), sd = .8), 6))

combs = c("10", "11", "20", "21", "22", "30", "31", "32", "33")

for (comb in combs) {
  assign(paste("ipol", comb, sep = "_"), interpolation(evalfun = evalfun, grid_list = get(paste("grid_list", comb, sep = "_")), cores = 80, comb = comb, tol = 1e-8, ratio = .9)$interpolant)
}
save(list = c(paste("ipol", combs, sep = "_")), file = "interpolant.rda", compress = "xz")
#usethis::use_data(ipol_10, ipol_11, ipol_20, ipol_21, ipol_22, ipol_30, ipol_31, ipol_32, ipol_33, internal = TRUE, compress = "xz")
