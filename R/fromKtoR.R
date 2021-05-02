##### This is going to be approximate version of fromKtoR which is much faster using multilinear interpolation. (07-23-2020)
##### multilinear interpolation method is from ipol in chebpol package.

# K: Kendall's tau matrix.
# zratio: a column vector of zero proportion values.
fromKtoR <- function(K, zratio = NULL, type = "trunc", method = "approx", tol = 1e-3) {
  K <- as.matrix(K)
  d1 <- nrow(K)
  p <- ifelse(is.null(ncol(zratio)), 1, ncol(zratio))
  # If this is just 1 variable, then correlation is automatically 1
  if (d1 == 1){return(as.matrix(1))}
  if (type == "continuous") {
    hatR <- sin(pi/2 * K)
  } else { # if the type is either "trunc" or "binary"
    upperR <- c(upper.tri(K)) # length p^2 of true/false with true corresponding to upper.tri
    hatRupper <- rep(NA, sum(upperR)) # length p(p-1)/2
    Kupper <- c(K[upperR]) # upper triangle of K matrix
    zratio1mat <- zratio2mat <- NULL

    # check if there is any element that is outside of the safe boundary for interpolation.
    for (i in 1:p) {
      zratio1mat = cbind(zratio1mat, rep(zratio[ , i], d1)[upperR]) # length p(p-1)/2
      zratio2mat = cbind(zratio2mat, rep(zratio[ , i], each = d1)[upperR]) # length p(p-1)/2
    }
    cutoff <- cutoff(type1 = type, type2 = type, tau = abs(Kupper), zratio1 = zratio1mat, zratio2 = zratio2mat, method = method)
    # Interpolate only those elements that are inside
    hatRupper[which(!(cutoff))] <- r_sol(type1 = type, type2 = type, tau = Kupper[which(!(cutoff))], zratio1 = zratio1mat[which(!(cutoff)), ], zratio2 = zratio2mat[which(!(cutoff)), ], method = "ml")
    for(ind in which(cutoff)){
      hatRupper[ind] = r_sol(type1 = type, type2 = type, tau = Kupper[ind], zratio1 = zratio1mat[ind, ], zratio2 = zratio2mat[ind, ], method = "original", tol = tol)
    }
    # Get upperR into hatR
    hatR <- matrix(0, d1, d1)
    hatR[upperR] <- hatRupper
    hatR <- hatR + t(hatR)
    diag(hatR) <- 1
  }
  return(hatR)
}

# K12: Kendall's tau matrix.
# zratio1: a vector of zero proportion values for row variables. The length should match with nrow of K12.
# zratio2: a vector of zero proportion values for column variables. The length should match with ncol of K12.
fromKtoR_mixed <- function(K12, zratio1 = NULL, zratio2 = NULL, type1 = "trunc", type2 = "continuous", method = "approx", tol = 1e-3) {

  K12 <- as.matrix(K12)
  d1 <- nrow(K12)
  d2 <- ncol(K12)
  p1 <- ifelse(is.null(ncol(zratio1)), 1, ncol(zratio1))
  p2 <- ifelse(is.null(ncol(zratio2)), 1, ncol(zratio2))
  if (type1 == "continuous" & type2 == "continuous") {
    hatR <- sin(pi/2 * K12)
  } else {
    # if the case is either CT, TC, TT, BC, BB or TB.
    hatR <- matrix(NA, d1, d2)
    zratio1mat <- zratio2mat <- NULL
    # check if there is any element that is outside of the safe boundary for interpolation.
    for (i in 1:p1) {
      zratio1mat = cbind(zratio1mat, rep(zratio1[ , i], d2))
    }
    for (j in 1:p2) {
      zratio2mat = cbind(zratio2mat, rep(zratio2[ , j], each = d1))
    }
    cutoff <- cutoff(type1 = type1, type2 = type2, tau = abs(K12), zratio1 = zratio1mat, zratio2 = zratio2mat, method = method)
    # Interpolate only those elements that are inside
    hatR[which(!(cutoff))] <- r_sol(type1 = type1, type2 = type2, tau = c(K12[which(!(cutoff))]), zratio1 = zratio1mat[which(!(cutoff)), ], zratio2 = zratio2mat[which(!(cutoff)), ], method = "ml")
    for(ind in which(cutoff)){
      hatR[ind] = r_sol(type1 = type1, type2 = type2, tau = K12[ind], zratio1 = zratio1mat[ind, ], zratio2 = zratio2mat[ind, ], method = "original", tol = tol)
    }
    # done with for the pairs that are outside of the safe boundary for multi-linear interpolation.
  }
  return(hatR)
}
