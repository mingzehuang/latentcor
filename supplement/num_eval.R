## The following script generate interpolant on my HPRC. Users should modify source directory and cores to replicate it.

library(devtools)
library(usethis)
library(stats)
library(fMultivar)
library(mnormt)
library(utils)
library(parallel)
library(future)
library(MASS)
library(foreach)
library(Matrix)
library(pcaPP)
library(doRNG)
library(doFuture)
library(heatmaply)
library(ggplot2)
library(plotly)
library(graphics)
library(microbenchmark)
library(latentcor)

#load("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/sysdata.rda")
#source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/internal.R")
#source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/GenData.R")
#source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/estR.R")
#source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/evaluation.R")

grid_list = list(LatentR = seq(-0.9, 0.9, by = 0.1), TruncRate = seq(0.1, 0.9, by = 0.1))

## For BC case
genfun = function (grid_input, ...) {
  out = gen_data(rhos = grid_input[1], XP = list(grid_input[2], NA), ...)$X
  return(out)
}
estfun_1 = function(X, ...) {
  out = latentcor(X = X, method = "original", ...)$R[1, 2]
  return(out)
}
estfun_2 = function(X, ...) {
  out = latentcor(X = X, method = "approx", ...)$R[1, 2]
  return(out)
}
evaluation_BC = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE, cores = 72, types = c("bin", "con"))

## For BB case
genfun = function (grid_input, ...) {
  out = gen_data(rhos = grid_input[1], XP = list(grid_input[2], .5), ...)$X
  return(out)
}
estfun_1 = function(X, ...) {
  out = latentcor(X = X, method = "original", ...)$R[1, 2]
  return(out)
}
estfun_2 = function(X, ...) {
  out = latentcor(X = X, method = "approx", ...)$R[1, 2]
  return(out)
}
evaluation_BB = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE, cores = 72, types = c("bin", "bin"))

## For TC case
genfun = function (grid_input, ...) {
  out = gen_data(rhos = grid_input[1], XP = list(grid_input[2], NA), ...)$X
  return(out)
}
estfun_1 = function(X, ...) {
  out = latentcor(X = X, method = "original", ...)$R[1, 2]
  return(out)
}
estfun_2 = function(X, ...) {
  out = latentcor(X = X, method = "approx", ...)$R[1, 2]
  return(out)
}
evaluation_TC = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE, cores = 72, types = c("tru", "con"))

## For TB case
genfun = function (grid_input, ...) {
  out = gen_data(rhos = grid_input[1], XP = list(grid_input[2], 0.5), ...)$X
  return(out)
}
estfun_1 = function(X, ...) {
  out = latentcor(X = X, method = "original", ...)$R[1, 2]
  return(out)
}
estfun_2 = function(X, ...) {
  out = latentcor(X = X, method = "approx", ...)$R[1, 2]
  return(out)
}
evaluation_TB = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE, cores = 72, types = c("tru", "bin"))

## For TT case
genfun = function (grid_input, ...) {
  out = gen_data(rhos = grid_input[1], XP = list(grid_input[2], 0.5), ...)$X
  return(out)
}
estfun_1 = function(X, ...) {
  out = latentcor(X = X, method = "original", ...)$R[1, 2]
  return(out)
}
estfun_2 = function(X, ...) {
  out = latentcor(X = X, method = "approx", ...)$R[1, 2]
  return(out)
}
evaluation_TT = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE, cores = 72, types = c("tru", "tru"))

## For NC case
genfun = function (grid_input, ...) {
  out = gen_data(rhos = grid_input[1], XP = list(c(grid_input[2], (1 - grid_input[2]) / 2), NA), ...)$X
  return(out)
}
estfun_1 = function(X, ...) {
  out = latentcor(X = X, method = "original", ...)$R[1, 2]
  return(out)
}
estfun_2 = function(X, ...) {
  out = latentcor(X = X, method = "approx", ...)$R[1, 2]
  return(out)
}
evaluation_NC = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE, cores = 72, types = c("ter", "con"))

## For NB case
genfun = function (grid_input, ...) {
  out = gen_data(rhos = grid_input[1], XP = list(c(grid_input[2], (1 - grid_input[2]) / 2), .5), ...)$X
  return(out)
}
estfun_1 = function(X, ...) {
  out = latentcor(X = X, method = "original", ...)$R[1, 2]
  return(out)
}
estfun_2 = function(X, ...) {
  out = latentcor(X = X, method = "approx", ...)$R[1, 2]
  return(out)
}
evaluation_NB = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE, cores = 72, types = c("ter", "bin"))

## For NT case
genfun = function (grid_input, ...) {
  out = gen_data(rhos = grid_input[1], XP = list(c(grid_input[2], (1 - grid_input[2]) / 2), .5), ...)$X
  return(out)
}
estfun_1 = function(X, ...) {
  out = latentcor(X = X, method = "original", ...)$R[1, 2]
  return(out)
}
estfun_2 = function(X, ...) {
  out = latentcor(X = X, method = "approx", ...)$R[1, 2]
  return(out)
}
evaluation_NT = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE, cores = 72, types = c("ter", "tru"))

## For NN case
genfun = function (grid_input, ...) {
  out = gen_data(rhos = grid_input[1], XP = list(c(grid_input[2], (1 - grid_input[2]) / 2), c(grid_input[2] / 2, (1 - grid_input[2]) / 4)), ...)$X
  return(out)
}
estfun_1 = function(X, ...) {
  out = latentcor(X = X, method = "original", ...)$R[1, 2]
  return(out)
}
estfun_2 = function(X, ...) {
  out = latentcor(X = X, method = "approx", ...)$R[1, 2]
  return(out)
}
evaluation_NN = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE, cores = 72, types = c("ter", "ter"))
save(evaluation_BC, evaluation_BB, evaluation_TC, evaluation_TB, evaluation_TT, evaluation_NC, evaluation_NB, evaluation_NT, evaluation_NN, file = "all_evaluation.rda", compress = "xz")

