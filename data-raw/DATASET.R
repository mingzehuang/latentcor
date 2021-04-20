## code to prepare `DATASET` dataset goes here

# usethis::use_data("DATASET")

# grid is revised with coarser grid on August 20, 2020.

############################################################################################
# For multilinear interpolation approximation for bridge Inverse
############################################################################################
require(foreach)
require(doParallel)
library(foreach)
library(doParallel)

# For BC case
############################################################################################
# load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/bc_0817.Rda")
#
# # grid values that used to create precomputed values
# # d1 <- seq(0.01, 0.99, length.out = 50)
# # tau1 <- c(seq(-0.5, -0.1, by = 0.007), seq(-0.095, -0.001, by = 0.005))
# # tau <- c(tau1, 0, rev(-tau1))
#
# # create input values for ipol
# BCvalue <- matrix(unlist(gridBCinv), ncol = length(d1), byrow = FALSE)

# # create grid input for ipol
# BCipolgrid <- list(tau, d1)
#
# # interpolation
# BCipol <- chebpol::ipol(BCvalue, grid = BCipolgrid, method = "multilin")

# For BC Case
# grid values that used to create precomputed values
tau1_grid <- c(seq(-0.5, -0.1, by = 0.007), seq(-0.095, -0.001, by = 0.005))
tau_grid <- c(tau1_grid, 0, rev(-tau1_grid))
d1_grid <- seq(0.02, 0.98, by = .02)
l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid)
BCvalue <- matrix(NA, l_tau_grid, l_d1_grid)

# # Single core single thread version
# for (i in 1:l_tau_grid) {
#   for (j in 1:l_d1_grid) {
#       f1 <- function(r)(bridgeF_bc(r, zratio1 = d1_grid[j]) - tau_grid[i])^2
#       op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
#       if(op == 100) {
#         warning("Optimize returned error one of the pairwise correlations, returning NA")
#         BCvalue[i, j] <- NA
#       } else {
#         BCvalue[i, j] <- unlist(op)
#       }
#   }
# }

# Single node multiple cores version
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tau_grid, .combine = rbind) %:%
  foreach (j = 1:l_d1_grid, .combine = c) %dopar% {
    f1 <- function(r)(bridgeF_bc(r, zratio1 = d1_grid[j]) - tau_grid[i])^2
    op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
    if(op == 100) {
      warning("Optimize returned error one of the pairwise correlations, returning NA")
      value_list <- NA
    } else {
      value_list <- unlist(op)
    }
  }
stopCluster(cl)

BCvalue = value_list

# create grid input for ipol
BCipolgrid <- list(tau_grid, d1_grid)
# interpolation.
BCipol <- chebpol::ipol(BCvalue, grid = BCipolgrid, method = "multilin")
save(BCipol, file = "BC_grid.rda", compress = "xz")
############################################################################################

# For BB case
############################################################################################
# load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/bb_0817.Rda")
#
# # grid values that used to create precomputed values
# # d1 <- d2 <- seq(0.01, 0.99, length.out = 50)
# # tau1 <- c(seq(-0.5, -0.1, by = 0.007), seq(-0.095, -0.001, by = 0.005))
# # tau <- c(tau1, 0, rev(-tau1))
#
# BBvalue <- array(NA, dim = c(length(tau), length(d1), length(d2)))
# for (i in 1:length(d1)){
#   for ( j in 1:length(d2)){
#     for ( k in 1:length(tau)){
#       BBvalue[k, i, j] <- gridBBinv[[length(d2)*(i - 1) + j]][k]
#     }
#   }
# }
#
# # create grid input for ipol
# BBipolgrid <- list(tau, d1, d2)
#
# # interpolation.
# BBipol <- chebpol::ipol(BBvalue, grid = BBipolgrid, method = "multilin")

# For BB Case
# grid values that used to create precomputed values
tau1_grid <- c(seq(-0.5, -0.1, by = 0.007), seq(-0.095, -0.001, by = 0.005))
tau_grid <- c(tau1_grid, 0, rev(-tau1_grid))
d1_grid <- d2_grid <- seq(0.02, 0.98, by = 0.02)
l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid); l_d2_grid <- length(d2_grid)
BBvalue <- array(NA, c(l_tau_grid, l_d1_grid, l_d2_grid))

# # Single core single thread version
# for (i in 1:l_tau_grid) {
#   for (j in 1:l_d1_grid) {
#     for (k in 1:l_d2_grid) {
#       f1 <- function(r)(bridgeF_bb(r, zratio1 = d1_grid[j], zratio2 = d2_grid[k]) - tau_grid[i])^2
#       op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
#       if(op == 100) {
#         warning("Optimize returned error one of the pairwise correlations, returning NA")
#         BBvalue[i, j, k] <- NA
#       } else {
#         BBvalue[i, j, k] <- unlist(op)
#       }
#     }
#   }
# }

# Parallel version
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tau_grid) %:%
  foreach (j = 1:l_d1_grid, .combine = rbind) %dopar% {
    value <- rep(NA, l_d2_grid)
    for (k in 1:l_d2_grid) {
      f1 <- function(r)(bridgeF_bb(r, zratio1 = d1_grid[j], zratio2 = d2_grid[k]) - tau_grid[i])^2
      op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
      if(op == 100) {
        warning("Optimize returned error one of the pairwise correlations, returning NA")
        value[k] <- NA
      } else {
        value[k] <- unlist(op)
      }
    }
    value_list <- value
  }
stopCluster(cl)

for (i in 1:l_tau_grid) {
  BBvalue[i, , ] = value_list[[i]]
}

# create grid input for ipol
BBipolgrid <- list(tau_grid, d1_grid, d2_grid)
# interpolation.
BBipol <- chebpol::ipol(BBvalue, grid = BBipolgrid, method = "multilin")
save(BBipol, file = "BB_grid.rda")

############################################################################################
# For TC case
############################################################################################

# load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/tc_0804.Rda")
#
# # grid values that used to create precomputed values.
# # d1 <- log10(seq(1, 10^0.99, length = 50))
# # tau <- seq(-0.99, 0.99, by = 0.01) # "by" increased from 0.005 to 0.01.
#
# # create computed values (in matrix) and grid (in list) for ipol function.
# value <- matrix(unlist(gridTCinv), ncol = length(d1), byrow = FALSE)
# grid <- list(tau, d1) # the length of list should be the same as the kinds of inputs.
#
# interp_multilin <- chebpol::ipol(value, grid = grid, method = "multilin")
#
# # create input values for ipol
# TCvalue <- matrix(unlist(gridTCinv), ncol = length(d1), byrow = FALSE)
#
# # create grid input for ipol
# TCipolgrid <- list(tau, d1)
#
# # interpolation
# TCipol <- chebpol::ipol(TCvalue, grid = TCipolgrid, method = "multilin")

# For TC Case
# grid values that used to create precomputed values
tau_grid <- seq(-0.98, 0.98, by = 0.02) # "by" increased from 0.005 to 0.01.
d1_grid <- log10(seq(10^0.02, 10^0.98, length = 50))
l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid)
TCvalue <- matrix(NA, l_tau_grid, l_d1_grid)

# # Single core single thread version
# for (i in 1:l_tau_grid) {
#   for (j in 1:l_d1_grid) {
#         f1 <- function(r)(bridgeF_tc(r, zratio1 = d1_grid[j]) - tau_grid[i])^2
#         op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
#         if(op == 100) {
#           warning("Optimize returned error one of the pairwise correlations, returning NA")
#           TCvalue[i, j] <- NA
#         } else {
#           TCvalue[i, j] <- unlist(op)
#         }
#   }
# }

# Parallel version
cl <- makePSOCKcluster(detectCores()) # Use makeForkCluster() on Linux-based system.
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tau_grid, .combine = rbind) %:%
    foreach (j = 1:l_d1_grid, .combine = c) %dopar% {
      f1 <- function(r)(bridgeF_tc(r, zratio1 = d1_grid[j]) - tau_grid[i])^2
      op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
      if(op == 100) {
        warning("Optimize returned error one of the pairwise correlations, returning NA")
        value_list <- NA
      } else {
        value_list <- unlist(op)
      }
    }
stopCluster(cl)

TCvalue = value_list


# create grid input for ipol
TCipolgrid <- list(tau_grid, d1_grid)
# interpolation.
TCipol <- chebpol::ipol(TCvalue, grid = TCipolgrid, method = "multilin")
save(TCipol, file = "TC_grid.rda")

# For TB case
############################################################################################
# load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/tb_0817.Rda")
#
# # grid values that used to create precomputed values
# # d1 <- log10(seq(1, 10^0.99, length = 50))
# # d2 <- seq(0.01, 0.99, length.out = 50)
# # tau1 <- c(seq(-0.5, -0.1, by = 0.007), seq(-0.095, -0.001, by = 0.005))
# # tau <- c(tau1, 0, rev(-tau1))
#
# TBvalue <- array(NA, dim = c(length(tau), length(d1), length(d2)))
# for (i in 1:length(d1)){
#   for ( j in 1:length(d2)){
#     for ( k in 1:length(tau)){
#       TBvalue[k, i, j] <- gridTBinv[[length(d2)*(i - 1) + j]][k]
#     }
#   }
# }
#
# # create grid input for ipol
# TBipolgrid <- list(tau, d1, d2)
#
# # interpolation.
# TBipol <- chebpol::ipol(TBvalue, grid = TBipolgrid, method = "multilin")

# For TB Case
# grid values that used to create precomputed values
tau1_grid <- c(seq(-0.5, -0.1, by = 0.007), seq(-0.095, -0.001, by = 0.005))
tau_grid <- c(tau1_grid, 0, rev(-tau1_grid))
d1_grid <- log10(seq(10^0.02, 10^0.99, length = 50))
d2_grid <- seq(0.02, 0.98, by = 0.02)
l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid); l_d2_grid <- length(d2_grid)
TBvalue <- array(NA, c(l_tau_grid, l_d1_grid, l_d2_grid))

# # Single core single thread version
# for (i in 1:l_tau_grid) {
#   for (j in 1:l_d1_grid) {
#     for (k in 1:l_d2_grid) {
#       f1 <- function(r)(bridgeF_tb(r, zratio1 = d1_grid[j], zratio2 = d2_grid[k]) - tau_grid[i])^2
#       op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
#       if(op == 100) {
#         warning("Optimize returned error one of the pairwise correlations, returning NA")
#         TBvalue[i, j, k] <- NA
#       } else {
#         TBvalue[i, j, k] <- unlist(op)
#       }
#     }
#   }
# }

# Parallel version
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tau_grid) %:%
  foreach (j = 1:l_d1_grid, .combine = rbind) %dopar% {
    value = rep(NA, l_d2_grid)
    for (k in 1:l_d2_grid) {
      f1 <- function(r)(bridgeF_tb(r, zratio1 = d1_grid[j], zratio2 = d2_grid[k]) - tau_grid[i])^2
      op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
      if(op == 100) {
        warning("Optimize returned error one of the pairwise correlations, returning NA")
        value[k] <- NA
      } else {
        value[k] <- unlist(op)
      }
    }
    value_list <- value
  }
stopCluster(cl)

for (i in 1:l_tau_grid) {
  TBvalue[i, , ] = value_list[[i]]
}

# create grid input for ipol
TBipolgrid <- list(tau_grid, d1_grid, d2_grid)
# interpolation.
TBipol <- chebpol::ipol(TBvalue, grid = TBipolgrid, method = "multilin")
save(TBipol, file = "TB_grid.rda")

############################################################################################

############################################################################################
# For TT case
############################################################################################

# load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/tt_0804.Rda")
#
# # grid values that used to create precomputed values.
# # d2 <- log10(seq(1, 10^0.99, length = 50))
# # tau <- seq(-0.99, 0.99, by = 0.01) # "by" increased from 0.005 to 0.01.
#
# TTvalue <- array(NA, dim = c(length(tau), length(d1), length(d2)))
# for (i in 1:length(d1)){
#   for ( j in 1:length(d2)){
#     for ( k in 1:length(tau)){
#       TTvalue[k, i, j] <- gridTTinv[[length(d2)*(i - 1) + j]][k]
#     }
#   }
# }
#
# # create grid input for ipol
# TTipolgrid <- list(tau, d1, d2)
#
# # interpolation.
# TTipol <- chebpol::ipol(TTvalue, grid = TTipolgrid, method = "multilin")

# For TT Case
# grid values that used to create precomputed values
tau_grid <- seq(-0.98, 0.98, by = 0.02) # "by" increased from 0.005 to 0.01.
d1_grid <- d2_grid <- log10(seq(10^0.02, 10^0.98, length = 50))
l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid); l_d2_grid <- length(d2_grid)
TTvalue <- array(NA, c(l_tau_grid, l_d1_grid, l_d2_grid))

# # Single core single thread version
# for (i in 1:l_tau_grid) {
#   for (j in 1:l_d1_grid) {
#     for (k in 1:l_d2_grid) {
#           f1 <- function(r)(bridgeF_tt(r, zratio1 = d1_grid[j], zratio2 = d2_grid[k]) - tau_grid[i])^2
#           op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
#           if(op == 100) {
#             warning("Optimize returned error one of the pairwise correlations, returning NA")
#             TTvalue[i, j, k] <- NA
#           } else {
#             TTvalue[i, j, k] <- unlist(op)
#           }
#     }
#   }
# }

# Parallel version
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tau_grid) %:%
    foreach (j = 1:l_d1_grid, .combine = rbind) %dopar% {
      value = rep(NA, l_d2_grid)
      for (k in 1:l_d2_grid) {
        f1 <- function(r)(bridgeF_tt(r, zratio1 = d1_grid[j], zratio2 = d2_grid[k]) - tau_grid[i])^2
        op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
        if(op == 100) {
          warning("Optimize returned error one of the pairwise correlations, returning NA")
          value[k] <- NA
        } else {
          value[k] <- unlist(op)
        }
      }
      value_list <- value
    }
stopCluster(cl)

for (i in 1:l_tau_grid) {
      TTvalue[i, , ] = value_list[[i]]
}

# create grid input for ipol
TTipolgrid <- list(tau_grid, d1_grid, d2_grid)
# interpolation.
TTipol <- chebpol::ipol(TTvalue, grid = TTipolgrid, method = "multilin")
save(TTipol, file = "TT_grid.rda")
############################################################################################

# For NC Case
# grid values that used to create precomputed values
tau_grid <- seq(-0.99, 0.99, by = 0.01) # "by" increased from 0.005 to 0.01.
d11_grid <- d12_grid <- seq(0.01, 0.99, length.out = 50)
l_tau_grid <- length(tau_grid); l_d11_grid <- length(d11_grid); l_d12_grid <- length(d12_grid)
NCvalue <- array(NA, c(l_tau_grid, l_d11_grid, l_d12_grid))

# # Single core single thread version
# for (i in 1:l_tau_grid) {
#   for (j in 1:l_d11_grid) {
#     for (k in j:l_d12_grid) {
#         f1 <- function(r)(bridgeF_nc(r, zratio1 = c(d11_grid[j], d12_grid[k])) - tau_grid[i])^2
#         op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
#         if(op == 100) {
#           warning("Optimize returned error one of the pairwise correlations, returning NA")
#           NCvalue[i, j, k] <- NA
#         } else {
#           NCvalue[i, j, k] <- unlist(op)
#         }
#     }
#   }
# }

# Parallel version
cl <- makePSOCKcluster(detectCores(logical=T))
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tau_grid) %:%
  foreach (j = 1:l_d11_grid, .combine = rbind) %dopar% {
    value <- rep(NA, l_d12_grid)
    for (k in j:l_d12_grid) {
      f1 <- function(r)(bridgeF_nc(r, zratio1 = c(d11_grid[j], d12_grid[k])) - tau_grid[i])^2
      op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
      if(op == 100) {
        warning("Optimize returned error one of the pairwise correlations, returning NA")
        value[k] <- NA
      } else {
        value[k] <- unlist(op)
      }
    }
    value_list <- value
  }
stopCluster(cl)

for (i in 1:l_tau_grid) {
  NCvalue[i, , ] = matrix(as.integer(10^7 * value_list[[i]]), l_d11_grid, l_d12_grid)
}

# create grid input for ipol
NCipolgrid <- list(tau_grid, d11_grid, d12_grid)
# interpolation.
NCipol <- chebpol::ipol(NCvalue, grid = NCipolgrid, method = "multilin")
save(NCipol, file = "NC_grid.rda", compress = "xz")


# For NB Case
# grid values that used to create precomputed values
tau_grid <- seq(-0.98, 0.98, by = 0.02) # "by" increased from 0.005 to 0.01.
d_combn <- matrix(seq(0.02, 0.98, by = 0.02)[combn(49, 2)], nrow = 2)
d1_grid <- 1:ncol(d_combn)
d2_grid <- seq(0.02, 0.98, by = 0.02)
l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid)
l_d2_grid <- length(d2_grid)
NBvalue <- array(NA, c(l_tau_grid, l_d1_grid, l_d2_grid))

# # Single core single thread version
# for (i in 1:l_tau_grid) {
#   for (j in 1:l_d11_grid) {
#     for (k in j:l_d12_grid) {
#       for (l in 1:l_d2_grid) {
#           f1 <- function(r)(bridgeF_nb(r, zratio1 = c(d11_grid[j], d12_grid[k]), zratio2 = d2_grid[l]) - tau_grid[i])^2
#           op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
#           if(op == 100) {
#             warning("Optimize returned error one of the pairwise correlations, returning NA")
#             NBvalue[i, j, k, l] <- NA
#           } else {
#             NBvalue[i, j, k, l] <- unlist(op)
#           }
#       }
#     }
#   }
# }

# Parallel version
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tau_grid) %:%
  foreach (j = 1:l_d1_grid, .combine = rbind) %dopar% {
    value = rep(NA, l_d2_grid)
    for (k in 1:l_d2_grid) {
      f1 <- function(r)(bridgeF_nb(r, zratio1 = c(d_combn[1, j], d_combn[2, j]), zratio2 = d2_grid[k]) - tau_grid[i])^2
      op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
      if(op == 100) {
        warning("Optimize returned error one of the pairwise correlations, returning NA")
        value[k] <- NA
      } else {
        value[k] <- unlist(op)
      }
    }
    value_list <- value
  }
stopCluster(cl)

for (i in 1:l_tau_grid) {
  NBvalue[i, , ] = matrix(as.integer(10^7 * value_list[[i]]), l_d1_grid, l_d2_grid)
}

# create grid input for ipol
NBipolgrid <- list(tau_grid, d1_grid, d2_grid)
# interpolation.
NBipol <- chebpol::ipol(NBvalue, grid = NBipolgrid, method = "multilin")
save(NBipol, file = "NB_grid.rda")


# For NN Case
# grid values that used to create precomputed values
tau_grid <- seq(-0.98, 0.98, by = 0.02) # "by" increased from 0.005 to 0.01.
d_combn <- matrix(seq(0.02, 0.98, by = 0.02)[combn(49, 2)], nrow = 2)
d1_grid <- d2_grid <- 1:ncol(d_combn);
l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid)
l_d2_grid <- length(d2_grid)
NNvalue <- array(NA, c(l_tau_grid, l_d1_grid, l_d2_grid))

# # Single core single thread version
# for (i in 1:l_tau_grid) {
#   for (j in 1:l_d11_grid) {
#     for (k in j:l_d12_grid) {
#       for (l in 1:l_d21_grid) {
#         for (m in l:l_d22_grid) {
#           f1 <- function(r)(bridgeF_nn(r, zratio1 = c(d11_grid[j], d12_grid[k]), zratio2 = c(d21_grid[l], d22_grid[m])) - tau_grid[i])^2
#           op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
#           if(op == 100) {
#             warning("Optimize returned error one of the pairwise correlations, returning NA")
#             NNvalue[i, j, k, l, m] <- NA
#           } else {
#             NNvalue[i, j, k, l, m] <- unlist(op)
#           }
#         }
#       }
#     }
#   }
# }

# Parallel version
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tau_grid) %:%
    foreach (j = 1:l_d1_grid, .combine = rbind) %dopar% {
      value = rep(NA, l_d2_grid)
      for (k in j:l_d2_grid) {
        f1 <- function(r)(bridgeF_nn(r, zratio1 = c(d_combn[1, j], d_combn[2, j]),
                        zratio2 = c(d_combn[1, k], d_combn[2, k])) - tau_grid[i])^2
        op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = 1e-3)[1], error = function(e) 100)
        if(op == 100) {
          warning("Optimize returned error one of the pairwise correlations, returning NA")
          value[k] <- NA
        } else {
          value[k] <- unlist(op)
        }
      }
      value_list <- value
    }
stopCluster(cl)

for (i in 1:l_tau_grid) {
  NNvalue[i, , ] = matrix(as.integer(10^7 * value_list[[i]]), l_d1_grid, l_d2_grid)
}
# create grid input for ipol
NNipolgrid <- list(tau_grid, d1_grid, d2_grid)
# interpolation.
NNipol <- chebpol::ipol(NNvalue, grid = NNipolgrid, method = "multilin")
save(NNipol, file = "NN_grid.rda")

usethis::use_data(BCipol, BBipol, TCipol, TBipol, TTipol, NCipol, NBipol, NNipol, internal = TRUE, overwrite = TRUE, compress = "xz")

library(MASS)
library(microbenchmark)
# setup for 100 replication.
nrep <- 100

# sample size
n <- 100

# will test 9 latent r and 11 zero proportion values.
latentRseq <- seq(0.05, 0.91, length.out = 9)
zratioseq <- c(0.04, 0.16, 0.28, 0.36, 0.44, 0.5, 0.56, 0.64, 0.72, 0.84, 0.96)

##### check six cases of TC, TT, BC, BB, TB
for (cases in 1:5){

  if(cases == 1){
    type1 <- "trunc"; type2 <- "continuous"
    typesh <- "TC"
  } else if(cases == 2){
    type1 <- "trunc"; type2 <- "trunc"
    typesh <- "TT"
  } else if(cases == 3){
    type1 <- "binary"; type2 <- "continuous"
    typesh <- "BC"
  } else if(cases == 4){
    type1 <- "binary"; type2 <- "binary"
    typesh <- "BB"
  } else if(cases == 5){
    type1 <- "trunc"; type2 <- "binary"
    typesh <- "TB"
  }

  # the computation results will be saved in data.frame format
  df_comptime <- df_accuracy <- NULL


  for (trueR in latentRseq){

    for (zrate in zratioseq){
      # initialize for every combination
      Kcor_org <- Kcor_ml <- Kcor_mlbd <- rep(NA, nrep)
      time_org <- time_ml <- time_mlbd <- rep(NA, nrep)
      time_all <- matrix(NA, nrow = nrep, ncol = 3)

      ptm <- proc.time()
      set.seed(123)

      for(i in 1:nrep){
        # generate bivariate normal
        z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, trueR, trueR, 1), nrow=2))
        if(cases == 1){ ## TC case
          # shifting to control the truncation levels
          z1shift <- quantile(z[, 1], zrate)
          z1 <- z[, 1] - z1shift
          z2 <- z[, 2] - z1shift # shift the same amount of z1.
          # truncate the first variable.
          u1 <- ifelse(z1 > 0, z1, 0)
          u2 <- z2 # since this is continuous variable
        } else if (cases == 2){ ## TT case
          z1shift <- quantile(z[, 1], zrate) # shifting to control the truncation levels
          z2shift <- quantile(z[, 2], zrate/2) # shifting half of zrate to control the truncation rate of 2nd variable as half trucation level of the first variable.
          z1 <- z[, 1] - z1shift
          z2 <- z[, 2] - z2shift
          # truncation first and second variables.
          u1 <- ifelse(z1 > 0, z1, 0)
          u2 <- ifelse(z2 > 0, z2, 0)
        } else if (cases == 3){ ## BC case
          z1shift <- quantile(z[, 1], zrate) # shifting to control the truncation levels
          z1 <- z[, 1] - z1shift
          z2 <- z[, 2] - z1shift # shift the same amount of z1.
          # truncation first
          u1 <- ifelse(z1 > 0, 1, 0)
          u2 <- z2 # since this is continuous variable
        } else if (cases == 4){ ## BB case
          z1shift <- quantile(z[, 1], zrate) # shifting to control the truncation levels
          z2shift <- quantile(z[, 2], 0.5) # fixed the zero ratio of the second variable
          z1 <- z[, 1] - z1shift
          z2 <- z[, 2] - z2shift
          # truncation first
          u1 <- ifelse(z1 > 0, 1, 0)
          u2 <- ifelse(z2 > 0, 1, 0)
        } else if (cases == 5){ ## TB case
          z1shift <- quantile(z[, 1], zrate) # shifting to control the truncation levels
          z2shift <- quantile(z[, 2], 0.5) # fixed the zero ratio of the second variable
          z1 <- z[, 1] - z1shift
          z2 <- z[, 2] - z2shift
          # truncation first
          u1 <- ifelse(z1 > 0, z1, 0)
          u2 <- ifelse(z2 > 0, 1, 0)
        }
        # didn't apply any transformation.
        x1 <- u1
        x2 <- u2


        capture.output( # suppress the microbenchmark result.
          time_all[i, ] <- print(microbenchmark(
            Kcor_org[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "original", nu = 0, tol = 1e-6)$R12,
            # Kcor_ml[i] <- estimateR_mixed_mlonly(X1 = x1, X2 = x2, type1 = type1, type2 = type2, nu = 0)$R12,
            Kcor_ml[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, nu = 0)$R12,
            Kcor_mlbd[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "approx", nu = 0, tol = 1e-6)$R12,
            times = 10 # tried ten times
          ), unit = "us")[, 5] # to use fixed unit: "microseconds"
          # 5th column has median value
        )
      }

      # apply(time_all, 2, summary) # 3rd row gives median value.
      df_comptime <- rbind.data.frame(df_comptime, data.frame(LatentR = trueR, TruncRate = zrate, medianTime = apply(time_all, 2, summary)[3, ], method = c("org", "ipol", "ipol_UB")))

      # save two kinds of errors: maximum absolute error and mean absolute error.
      df_accuracy <- rbind.data.frame(df_accuracy, data.frame(LatentR = trueR, TruncRate = zrate, MeanAD = c(mean(abs(Kcor_org - Kcor_ml)), mean(abs(Kcor_org - Kcor_mlbd))), MaxAD = c(max(abs(Kcor_org - Kcor_ml)), max(abs(Kcor_org - Kcor_mlbd))), method = c("ipol", "ipol_UB")))

      cat(typesh, "case: trueR = ", trueR, "\t zrate =", zrate, "\t took ", (proc.time() - ptm)[3], " seconds.\n")
    }
  }


  save(df_comptime, df_accuracy, file = paste0("TwoSim_", typesh, "_rep10.rda"))
}


library(MASS)
library(microbenchmark)
library(foreach)
library(doParallel)
# setup for 100 replication.
nrep <- 100

# sample size
n <- 100

# will test 9 latent r and 11 zero proportion values.
latentRseq <- seq(0.05, 0.91, length.out = 9)
zratioseq <- c(0.04, 0.16, 0.28, 0.36, 0.44, 0.5, 0.56, 0.64, 0.72, 0.84, 0.96)

##### check BC
type1 <- "binary"; type2 <- "binary"
typesh <- "BB"
# the computation results will be saved in data.frame format
df_comptime <- df_accuracy <- NULL

cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)
BB_eval <-
  foreach (trueR = 1:length(latentRseq)) %:%
  foreach (zrate = 1:length(zratioseq)) %dopar% {
    # initialize for every combination
    Kcor_org <- Kcor_ml <- Kcor_mlbd <- rep(NA, nrep)
    time_org <- time_ml <- time_mlbd <- rep(NA, nrep)
    time_all <- matrix(NA, nrow = nrep, ncol = 3)
    set.seed(123)
    for(i in 1:nrep){
      # generate bivariate normal
      z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, latentRseq[trueR], latentRseq[trueR], 1), nrow=2))
      z1shift <- quantile(z[, 1], zratioseq[zrate]) # shifting to control the truncation levels
      z2shift <- quantile(z[, 2], 0.5) # fixed the zero ratio of the second variable
      z1 <- z[, 1] - z1shift
      z2 <- z[, 2] - z2shift
      # truncation first
      u1 <- ifelse(z1 > 0, 1, 0)
      u2 <- ifelse(z2 > 0, 1, 0)
      # didn't apply any transformation.
      x1 <- u1
      x2 <- u2
      capture.output( # suppress the microbenchmark result.
        time_all[i, ] <- print(microbenchmark::microbenchmark(
          Kcor_org[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "original", nu = 0, tol = 1e-6)$R12,
          Kcor_ml[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, nu = 0)$R12,
          Kcor_mlbd[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "approx", nu = 0, tol = 1e-6)$R12,
          times = 10 # tried ten times
        ), unit = "us")[, 5] # to use fixed unit: "microseconds"
        # 5th column has median value
      )
    }
    # apply(time_all, 2, summary) # 3rd row gives median value.
    BB_eval <- data.frame(LatentR = latentRseq[trueR], TruncRate = zratioseq[zrate], medianTime = apply(time_all, 2, summary)[3, ], MeanAD = c(0, mean(abs(Kcor_org - Kcor_ml)), mean(abs(Kcor_org - Kcor_mlbd))), MaxAD = c(0, max(abs(Kcor_org - Kcor_ml)), max(abs(Kcor_org - Kcor_mlbd))), method = c("org", "ipol", "ipol_UB"))
  }
save(BB_eval, file = "BB_eval")
save(BCipol, BBipol, TCipol, TBipol, TTipol, NCipol, NBipol, NNipol, file = "all_grid.rda",
     compress = "xz")
