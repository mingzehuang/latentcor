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

nrep = 1:100
types_comb = expand.grid(c("ter", "tru", "bin", "con"), c("ter", "tru", "bin", "con"))[upper.tri(matrix(1:16, 4, 4), diag = TRUE), ]
types_list = list()
for (i in 1:(nrow(types_comb) - 1)) {
  types_list[[i]] = types_comb[i, ]
}
latentRseq = seq(-0.9, 0.9, by = 0.1); zratioseq = seq(0.1, 0.9, by = 0.1); ratios = c(0, 0.5, 0.9, 0.99, 1)
indseq = expand.grid(nrep, types_list, latentRseq, zratioseq, ratios)
registerDoFuture()
plan(multicore, workers = 80)
value =
  foreach (ind = 1:nrow(indseq), .combine = c) %dopar% {
    types = unlist(indseq[ind, 2])
    if (types[1] == "ter"  & types[2] == "ter") {
      XP = list(c(indseq[ind, 4] / 2, indseq[ind, 4]), c(indseq[ind, 4] / 2, indseq[ind, 4]))
    } else if (types[1] == "ter" & (types[2] == "tru" | types[2] == "bin")) {
      XP = list(c(indseq[ind, 4] / 2, indseq[ind, 4]), 0.5)
    } else if (types[1] == "ter" & types[2] == "con") {
      XP = list(c(indseq[ind, 4] / 2, indseq[ind, 4]), NA)
    } else if ((types[1] == "tru" | types[1] == "bin") & (types[2] == "tru" | types[2] == "bin")) {
      XP = list(indseq[ind, 4], 0.5)
    } else if ((types[1] == "tru" | types[1] == "bin") & type[2] == "con") {
      XP = list(indseq[ind, 4], NA)
    }
    simdata = GenData(types = types, rhos = indseq[ind, 3], XP = XP)
    X = simdata$X
    time = median(microbenchmark::microbenchmark(estimate = estR(X = simdata$X, types = unlist(indseq[ind, 2]), ratio = indseq[ind, 5])$R[1, 2], times = 5)$time) / 10^6
    value = c(time, estimate)
  }
time_mat = array(value[ , 1], c(length(nrep), length(types_list), length(latentRseq), length(zratioseq), length(ratios)))
estimate_mat = array(value[ , 2], c(length(nrep), length(types_list), length(latentRseq), length(zratioseq), length(ratios)))

median_time = array(NA, c(length(types_list), length(latentRseq), length(zratioseq), length(ratios)))
for (j in 1:length(types_list)) {
  for (k in 1:length(latentRseq)) {
    for (l in 1:length(zratioseq)) {
      for (m in 1:length(ratios)) {
        median_time[j, k, l, m] = median(time_mat[ , j, k, l, m])
      }
    }
  }
}

meanAE = array(NA, c(length(types_list), length(latentRseq), length(zratioseq), length(ratios)))
for (j in 1:length(types_list)) {
  for (k in 1:length(latentRseq)) {
    for (l in 1:length(zratioseq)) {
      for (m in 1:length(ratios)) {
        meanAE[j, k, l, m] = mean(abs(estimate_mat[ , j, k, l, m] - latentRseq[k]))
      }
    }
  }
}

meanAAE = array(NA, c(length(types_list), length(latentRseq), length(zratioseq), length(ratios) - 1))
for (j in 1:length(types_list)) {
  for (k in 1:length(latentRseq)) {
    for (l in 1:length(zratioseq)) {
      for (m in 1:(length(ratios) - 1)) {
        meanAAE[j, k, l, m] = mean(abs(estimate_mat[ , j, k, l, m + 1] - estimate_mat[ , j, k, l, 1]))
      }
    }
  }
}
