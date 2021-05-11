##### This is going to be approximate version of fromKtoR which is much faster using multilinear interpolation. (07-23-2020)
##### multilinear interpolation method is from ipol in chebpol package.

# K: Kendall's tau matrix.
# zratio: a column vector of zero proportion values.
fromKtoR <- function(K, zratio, type, method, tol, ratio) {
  K = as.matrix(K)
  # If this is just 1 variable, then correlation is automatically 1
  if (length(K) == 1){return(as.matrix(1))}
  else {
    d1 <- nrow(K)
    zratio1 = apply(zratio, 2, function(x){rep(x, d1)})
    zratio2 = apply(zratio, 2, function(x){rep(x, each = d1)})
    upperR <- c(upper.tri(K)) # length p^2 of true/false with true corresponding to upper.tri
    hatRupper = R_sol(type1 = type, type2 = type, tau = c(K[upperR]), zratio1 = matrix(zratio1[upperR, ], nrow = sum(upperR)), zratio2 = matrix(zratio2[upperR, ], nrow = sum(upperR)), method = method, tol = tol, ratio = ratio)
    # Get upperR into hatR
    hatR <- matrix(0, d1, d1)
    hatR[upperR] <- hatRupper
    hatR <- hatR + t(hatR)
    diag(hatR) <- 1
    return(hatR)
  }
}

# K12: Kendall's tau matrix.
# zratio1: a vector of zero proportion values for row variables. The length should match with nrow of K12.
# zratio2: a vector of zero proportion values for column variables. The length should match with ncol of K12.
fromKtoR_mixed <- function(K12, zratio1, zratio2, type1, type2, method, tol, ratio) {
  K12 = as.matrix(K12)
  d1 <- nrow(K12);  d2 <- ncol(K12)
  if (!(is.null(zratio1))) {zratio1 = apply(zratio1, 2, function(x){rep(x, d2)})}
  if (!(is.null(zratio2))) {zratio2 = apply(zratio2, 2, function(x){rep(x, each = d1)})}
  if ((type1 == "continuous" & type2 == "binary") | (type1 == "continuous" & type2 == "trunc") |
      (type1 == "binary" & type2 == "trunc") | (type1 == "continuous" & type2 == "ternary") |
      (type1 == "binary" & type2 == "ternary") | (type1 == "trunc" & type2 == "ternary")) {
    hatR = R_sol(type1 = type2, type2 = type1, tau = c(t(K12)), zratio1 = matrix(zratio2, nrow = length(K12)), zratio2 = matrix(zratio1, nrow = length(K12)), method = method, tol = tol, ratio = ratio)
    return(matrix(hatR, d1, d2, byrow = TRUE))
  } else {
    hatR = R_sol(type1 = type1, type2 = type2, tau = c(K12), zratio1 = matrix(zratio1, nrow = length(K12)), zratio2 = matrix(zratio2, nrow = length(K12)), method = method, tol = tol, ratio = ratio)
    return(matrix(hatR, d1, d2))
  }
}
