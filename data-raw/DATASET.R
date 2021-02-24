## code to prepare `DATASET` dataset goes here

# usethis::use_data("DATASET")

# grid is revised with coarser grid on August 20, 2020.

############################################################################################
# For multilinear interpolation approximation for bridge Inverse
############################################################################################


############################################################################################
# For TC case
############################################################################################

load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/tc_0804.Rda")

# grid values that used to create precomputed values.
# d1 <- log10(seq(1, 10^0.99, length = 50))
# tau <- seq(-0.99, 0.99, by = 0.01) # "by" increased from 0.005 to 0.01.

# create computed values (in matrix) and grid (in list) for ipol function.
value <- matrix(unlist(gridTCinv), ncol = length(d1), byrow = FALSE)
grid <- list(tau, d1) # the length of list should be the same as the kinds of inputs.

interp_multilin <- chebpol::ipol(value, grid = grid, method = "multilin")

# create input values for ipol
TCvalue <- matrix(unlist(gridTCinv), ncol = length(d1), byrow = FALSE)

# create grid input for ipol
TCipolgrid <- list(tau, d1)

# interpolation
TCipol <- chebpol::ipol(TCvalue, grid = TCipolgrid, method = "multilin")

############################################################################################
# For TT case
############################################################################################

load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/tt_0804.Rda")

# grid values that used to create precomputed values.
# d2 <- log10(seq(1, 10^0.99, length = 50))
# tau <- seq(-0.99, 0.99, by = 0.01) # "by" increased from 0.005 to 0.01.

TTvalue <- array(NA, dim = c(length(tau), length(d1), length(d2)))
for (i in 1:length(d1)){
  for ( j in 1:length(d2)){
    for ( k in 1:length(tau)){
      TTvalue[k, i, j] <- gridTTinv[[length(d2)*(i - 1) + j]][k]
    }
  }
}

# create grid input for ipol
TTipolgrid <- list(tau, d1, d2)

# interpolation.
TTipol <- chebpol::ipol(TTvalue, grid = TTipolgrid, method = "multilin")

############################################################################################
# For TB case
############################################################################################
load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/tb_0817.Rda")

# grid values that used to create precomputed values
# d1 <- log10(seq(1, 10^0.99, length = 50))
# d2 <- seq(0.01, 0.99, length.out = 50)
# tau1 <- c(seq(-0.5, -0.1, by = 0.007), seq(-0.095, -0.001, by = 0.005))
# tau <- c(tau1, 0, rev(-tau1))

TBvalue <- array(NA, dim = c(length(tau), length(d1), length(d2)))
for (i in 1:length(d1)){
  for ( j in 1:length(d2)){
    for ( k in 1:length(tau)){
      TBvalue[k, i, j] <- gridTBinv[[length(d2)*(i - 1) + j]][k]
    }
  }
}

# create grid input for ipol
TBipolgrid <- list(tau, d1, d2)

# interpolation.
TBipol <- chebpol::ipol(TBvalue, grid = TBipolgrid, method = "multilin")

############################################################################################
# For BC case
############################################################################################
load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/bc_0817.Rda")

# grid values that used to create precomputed values
# d1 <- seq(0.01, 0.99, length.out = 50)
# tau1 <- c(seq(-0.5, -0.1, by = 0.007), seq(-0.095, -0.001, by = 0.005))
# tau <- c(tau1, 0, rev(-tau1))

# create input values for ipol
BCvalue <- matrix(unlist(gridBCinv), ncol = length(d1), byrow = FALSE)

# create grid input for ipol
BCipolgrid <- list(tau, d1)

# interpolation
BCipol <- chebpol::ipol(BCvalue, grid = BCipolgrid, method = "multilin")

############################################################################################
# For BB case
############################################################################################
load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/bb_0817.Rda")

# grid values that used to create precomputed values
# d1 <- d2 <- seq(0.01, 0.99, length.out = 50)
# tau1 <- c(seq(-0.5, -0.1, by = 0.007), seq(-0.095, -0.001, by = 0.005))
# tau <- c(tau1, 0, rev(-tau1))

BBvalue <- array(NA, dim = c(length(tau), length(d1), length(d2)))
for (i in 1:length(d1)){
  for ( j in 1:length(d2)){
    for ( k in 1:length(tau)){
      BBvalue[k, i, j] <- gridBBinv[[length(d2)*(i - 1) + j]][k]
    }
  }
}

# create grid input for ipol
BBipolgrid <- list(tau, d1, d2)

# interpolation.
BBipol <- chebpol::ipol(BBvalue, grid = BBipolgrid, method = "multilin")


# For NN Case
# grid values that used to create precomputed values
tau_grid <- seq(-0.99, 0.99, by = 0.01) # "by" increased from 0.005 to 0.01.
d11_grid <- d12_grid <- d21_grid <- d22_grid <- seq(0.01, 0.99, length.out = 50)
NNvalue <- array(NA, c(length(tau_grid), length(d11_grid), length(d12_grid), length(d21_grid), length(d22_grid)))

for (i in 1:length(tau_grid)) {
  for (j in 1:length(d11_grid)) {
    for (k in j:length(d12_grid)) {
      for (l in 1:length(d21_grid)) {
        for (m in l:length(d22_grid)) {
          tau <- tau_grid[i]; d1 <- c(d11_grid[j], d12_grid[k]); d2 <- c(d21_grid[l], d22_grid[m])
          f1 <- function(r)(bridgeF_nn(r, zratio1 = d1, zratio2 = d2) - tau)^2
          op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
          if(op == 100) {
            warning("Optimize returned error one of the pairwise correlations, returning NA")
            NNvalue[i, j, k, l, m] <- NA
          } else {
            NNvalue[i, j, k, l, m] <- unlist(op)
          }
        }
      }
    }
  }
}
# create grid input for ipol
NNipolgrid <- list(tau_grid, d11_grid, d12_grid, d21_grid, d22_grid)
# interpolation.
NNipol <- chebpol::ipol(NNvalue, grid = NNipolgrid, method = "multilin")

# # For NB Case
# grid values that used to create precomputed values
tau_grid <- seq(-0.99, 0.99, by = 0.01) # "by" increased from 0.005 to 0.01.
d11_grid <- d12_grid <- d2_grid <- seq(0.01, 0.99, length.out = 50)
NBvalue <- array(NA, c(length(tau_grid), length(d11_grid), length(d12_grid), length(d2_grid)))

for (i in 1:length(tau_grid)) {
  for (j in 1:length(d11_grid)) {
    for (k in j:length(d12_grid)) {
      for (l in 1:length(d2_grid)) {
          tau <- tau_grid[i]; d1 <- c(d11_grid[j], d12_grid[k]); d2 <- d2_grid[l]
          f1 <- function(r)(bridgeF_nb(r, zratio1 = d1, zratio2 = d2) - tau)^2
          op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
          if(op == 100) {
            warning("Optimize returned error one of the pairwise correlations, returning NA")
            NBvalue[i, j, k, l] <- NA
          } else {
            NBvalue[i, j, k, l] <- unlist(op)
          }
      }
    }
  }
}
# create grid input for ipol
NBipolgrid <- list(tau_grid, d11_grid, d12_grid, d2_grid)
# interpolation.
NBipol <- chebpol::ipol(NBvalue, grid = NBipolgrid, method = "multilin")
#
# # For NC Case
# grid values that used to create precomputed values
tau_grid <- seq(-0.99, 0.99, by = 0.01) # "by" increased from 0.005 to 0.01.
d11_grid <- d12_grid <- seq(0.01, 0.99, length.out = 50)
NCvalue <- array(NA, c(length(tau_grid), length(d11_grid), length(d12_grid)))

for (i in 1:length(tau_grid)) {
  for (j in 1:length(d11_grid)) {
    for (k in j:length(d12_grid)) {
        tau <- tau_grid[i]; d1 <- c(d11_grid[j], d12_grid[k])
        f1 <- function(r)(bridgeF_nc(r, zratio1 = d1) - tau)^2
        op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
        if(op == 100) {
          warning("Optimize returned error one of the pairwise correlations, returning NA")
          NCvalue[i, j, k] <- NA
        } else {
          NCvalue[i, j, k] <- unlist(op)
        }
    }
  }
}
# create grid input for ipol
NCipolgrid <- list(tau_grid, d11_grid, d12_grid)
# interpolation.
NCipol <- chebpol::ipol(NCvalue, grid = NCipolgrid, method = "multilin")

usethis::use_data(TCipol, TTipol, TBipol, BCipol, BBipol, internal = TRUE, overwrite = TRUE, compress = "xz")

