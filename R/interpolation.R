#Library packages you need for internal functions
library(foreach)
library(parallel)
library(doFuture)
library(stats)
library(pcaPP)
library(fMultivar)
library(mnormt)
library(chebpol)
library(Matrix)
library(utils)
library(MASS)
source('~/latentcor_git/latentcor/R/internal.R')

value = function(evalfun, grid_list, cores = detectCores()) {
  grid_all = expand.grid(grid_list)
  registerDoFuture()
  plan(multicore, workers = cores)
  value_vector =
    foreach (j = 1:nrow(grid_all), .combine = c) %dopar% {
      grid_input = as.numeric(grid_all[j, ])
      value_list = evalfun(grid_input = grid_input)
    }
  d_grid = length(grid_list)
  dim_value = NULL
  for (i in 1:d_grid) {
    dim_value = c(dim_value, length(grid_list[[i]]))
  }
  return (array(as.integer(10^7 * value_vector), dim = dim_value))
}

evalfun_10 = function(grid_input){
  comb = "10"; zratio1 = grid_input[2]; zratio2 = NA
  K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
  r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
}

evalfun_11 = function(grid_input){
  comb = "11"; zratio1 = grid_input[2]; zratio2 = grid_input[3]
  K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
  r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
}

evalfun_20 = function(grid_input){
  comb = "20"; zratio1 = grid_input[2]; zratio2 = NA
  K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
  r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
}

evalfun_21 = function(grid_input){
  comb = "21"; zratio1 = grid_input[2]; zratio2 = grid_input[3]
  K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
  r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
}

evalfun_22 = function(grid_input){
  comb = "22"; zratio1 = grid_input[2]; zratio2 = grid_input[3]
  K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
  r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
}

evalfun_30 = function(grid_input){
  comb = "30"; zratio1 = c(grid_input[2] * grid_input[3], grid_input[3]); zratio2 = NA
  K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
  r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
}

evalfun_31 = function(grid_input){
  comb = "31"; zratio1 = c(grid_input[2] * grid_input[3], grid_input[3]); zratio2 = grid_input[4]
  K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
  r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
}

evalfun_32 = function(grid_input){
  comb = "32"; zratio1 = c(grid_input[2] * grid_input[3], grid_input[3]); zratio2 = grid_input[4]
  K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
  r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
}

evalfun_33 = function(grid_input){
  comb = "33"; zratio1 = c(grid_input[2] * grid_input[3], grid_input[3]); zratio2 = grid_input[4:5]
  K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
  r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
}

grid_list_10 = list(seq(-0.99, 0.99, by = 0.02), seq(0.02, 0.98, by = 0.02))
grid_list_11 = list(seq(-0.99, 0.99, by = 0.03), seq(0.02, 0.98, by = 0.03), seq(0.02, 0.98, by = 0.03))
grid_list_20 = list(seq(-0.99, 0.99, by = 0.02), seq(0.02, 0.98, by = 0.02))
grid_list_21 = list(seq(-0.99, 0.99, by = 0.03), seq(0.02, 0.98, by = 0.03), seq(0.02, 0.98, by = 0.03))
grid_list_22 = list(seq(-0.99, 0.99, by = 0.03), seq(0.02, 0.98, by = 0.03), seq(0.02, 0.98, by = 0.03))
grid_list_30 = list(seq(-0.99, 0.99, by = 0.03), seq(0.02, 0.98, by = 0.03), seq(0.02, 0.98, by = 0.03))
grid_list_31 = list(round(seq(-0.95, 0.95, by = 0.05)), round(seq(0.05, 0.95, by = 0.05)), round(seq(0.05, 0.95, by = 0.05)), round(seq(0.05, 0.95, by = 0.05)))
grid_list_32 = list(round(seq(-0.95, 0.95, by = 0.05)), round(seq(0.05, 0.95, by = 0.05)), round(seq(0.05, 0.95, by = 0.05)), round(seq(0.05, 0.95, by = 0.05)))
grid_list_33 = list(seq(-0.9, 0.9, by = 0.1), seq(0.1, 0.9, by = 0.1), seq(0.1, 0.9, by = 0.1), seq(0.1, 0.9, by = 0.1), seq(0.1, 0.9, by = 0.1))


combs = c("10", "11", "20", "21", "22", "30", "31", "32", "33")
ipols = NULL
for (comb in combs) {
  assign(paste("ipol", comb, sep = "_"), chebpol::ipol(value(evalfun = paste("evalfun", comb, sep = "_"), grid_list = paste("grid_list", comb, sep = "_")), grid = paste("grid_list", comb, sep = "_"), method = "multilin"))
}
save(list = c(paste("ipol", combs, sep = "_")), file = "interpolation.rda", compress = "xz")


