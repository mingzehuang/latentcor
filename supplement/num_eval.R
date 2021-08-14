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
