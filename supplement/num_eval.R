## The following script generate interpolant on my HPRC. Users should modify source directory and cores to replicate it.

library(MASS)
library(microbenchmark)
library(foreach)
library(Matrix)
library(chebpol)
library(pcaPP)
library(doRNG)
library(doFuture)

load("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/sysdata.rda")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/internal.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/GenData.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/estR.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/evaluation.R")


latentR = seq(-0.9, 0.9, by = 0.1); zratios = seq(0.1, 0.9, by = 0.1)
grid_list = list(latentR, zratios)

## For BC case
genfun = function (rho, zratio) {
  out = GenData(types = c("bin", "con"), rhos = rho, XP = list(zratio, NA))$X
  return(out)
}
estfun_1 = function(X) {
  out = estR(X = X, types = c("bin", "con"), method = "original")$R[1, 2]
  return(out)
}
estfun_2 = function(X) {
  out = estR(X = X, types = c("bin", "con"))$R[1, 2]
  return(out)
}
evaluation_BC = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE)

## For BB case
genfun = function (rho, zratio) {
  out = GenData(types = c("bin", "bin"), rhos = rho, XP = list(zratio, .5))$X
  return(out)
}
estfun_1 = function(X) {
  out = estR(X = X, types = c("bin", "bin"), method = "original")$R[1, 2]
  return(out)
}
estfun_2 = function(X) {
  out = estR(X = X, types = c("bin", "bin"))$R[1, 2]
  return(out)
}
evaluation_BB = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE)

## For TC case
genfun = function (rho, zratio) {
  out = GenData(types = c("tru", "con"), rhos = rho, XP = list(zratio, NA))$X
  return(out)
}
estfun_1 = function(X) {
  out = estR(X = X, types = c("tru", "con"), method = "original")$R[1, 2]
  return(out)
}
estfun_2 = function(X) {
  out = estR(X = X, types = c("tru", "con"))$R[1, 2]
  return(out)
}
evaluation_TC = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE)

## For TB case
genfun = function (rho, zratio) {
  out = GenData(types = c("tru", "bin"), rhos = rho, XP = list(zratio, 0.5))$X
  return(out)
}
estfun_1 = function(X) {
  out = estR(X = X, types = c("tru", "bin"), method = "original")$R[1, 2]
  return(out)
}
estfun_2 = function(X) {
  out = estR(X = X, types = c("tru", "bin"))$R[1, 2]
  return(out)
}
evaluation_TB = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE)

## For TT case
genfun = function (rho, zratio) {
  out = GenData(types = c("tru", "tru"), rhos = rho, XP = list(zratio, 0.5))$X
  return(out)
}
estfun_1 = function(X) {
  out = estR(X = X, types = c("tru", "tru"), method = "original")$R[1, 2]
  return(out)
}
estfun_2 = function(X) {
  out = estR(X = X, types = c("tru", "tru"))$R[1, 2]
  return(out)
}
evaluation_TT = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE)

## For NC case
genfun = function (rho, zratio) {
  out = GenData(types = c("ter", "con"), rhos = rho, XP = list(c(zratio, (1 - zratio) / 2), NA))$X
  return(out)
}
estfun_1 = function(X) {
  out = estR(X = X, types = c("ter", "con"), method = "original")$R[1, 2]
  return(out)
}
estfun_2 = function(X) {
  out = estR(X = X, types = c("ter", "con"))$R[1, 2]
  return(out)
}
evaluation_NC = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE)

## For NB case
genfun = function (rho, zratio) {
  out = GenData(types = c("ter", "bin"), rhos = rho, XP = list(c(zratio, (1 - zratio) / 2), .5))$X
  return(out)
}
estfun_1 = function(X) {
  out = estR(X = X, types = c("ter", "bin"), method = "original")$R[1, 2]
  return(out)
}
estfun_2 = function(X) {
  out = estR(X = X, types = c("ter", "bin"))$R[1, 2]
  return(out)
}
evaluation_NB = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE)

## For NT case
genfun = function (rho, zratio) {
  out = GenData(types = c("ter", "tru"), rhos = rho, XP = list(c(zratio, (1 - zratio) / 2), .5))$X
  return(out)
}
estfun_1 = function(X) {
  out = estR(X = X, types = c("ter", "tru"), method = "original")$R[1, 2]
  return(out)
}
estfun_2 = function(X) {
  out = estR(X = X, types = c("ter", "tru"))$R[1, 2]
  return(out)
}
evaluation_NT = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE)

## For NN case
genfun = function (rho, zratio) {
  out = GenData(types = c("ter", "ter"), rhos = rho, XP = list(c(zratio, (1 - zratio) / 2), c(zratio / 2, (1 - zratio) / 4)))$X
  return(out)
}
estfun_1 = function(X) {
  out = estR(X = X, types = c("ter", "ter"), method = "original")$R[1, 2]
  return(out)
}
estfun_2 = function(X) {
  out = estR(X = X, types = c("ter", "ter"))$R[1, 2]
  return(out)
}
evaluation_NN = evaluation(genfun = genfun, estfun_1 = estfun_1, estfun_2 = estfun_2, grid_list = grid_list, showplot = TRUE)
save(evaluation_BC, evaluation_BB, evaluation_TC, evaluation_TB, evaluation_TT, evaluation_NC, evaluation_NB, evaluation_NT, evaluation_NN, file = "all_evaluation.rda", compress = "xz")
