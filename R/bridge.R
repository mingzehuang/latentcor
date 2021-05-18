#'
#' @importFrom chebpol ipol
#'
NULL

zratio = function(X, type) {
  if (type == "continuous") {
    zratio = NULL
  } else if (type == "trunc"){
    zratio <- as.matrix(colMeans(X == 0))
    # checking data type
    if(sum(X < 0) > 0) {
      stop("The data of truncated type contains negative values.")
    }
    # checking proportion of zero values
    if(sum(zratio) == 0){
      message("The data does not contain zeros. Consider changing the type to \"continuous\".")
    }
    if (sum(zratio == 1) > 0){
      stop("There are variables in the data that have only zeros. Filter those     variables before continuing. \n")
    }
  } else if (type == "binary") {
    zratio <- as.matrix(colMeans(X == 0))
    # checking data type
    if(sum(!(X %in% c(0, 1))) > 0) {
      stop("The data is not \"binary\".")
    }
    if (sum(zratio == 1) > 0 | sum(zratio == 0) > 0){
      stop("There are binary variables in the data that have only zeros or only ones. Filter those variables before continuing. \n")
    }
  } else if (type == "ternary") {
    zratio <- cbind(colMeans(X == 0), 1 - colMeans(X == 2))
    # } else if (type == "dtrunc") {
    #   zratio <- cbind(colMeans(X == 0), 1 - colMeans(X == 1))
  }
  return(zratio)
}


R_sol <- function(type1, type2, tau, zratio1, zratio2, method, tol, ratio) {
  out <- rep(NA, length(tau))
  cutoff <- cutoff(type1 = type1, type2 = type2, tau = abs(tau), zratio1 = zratio1, zratio2 = zratio2, method = method, ratio = ratio)
  outside <- which(cutoff)
  inside <- which(!(cutoff))
  if (length(inside) > 0) {
    out[inside] <- r_ml(type1 = type1, type2 = type2, tau = tau[inside], zratio1 = matrix(zratio1[inside, ], nrow=length(inside)), zratio2 = matrix(zratio2[inside, ], nrow = length(inside)))
  }
  if (length(outside) > 0) {
    out[outside] = sapply(outside, function(x){r_sol(type1 = type1, type2 = type2, tau = tau[x], zratio1 = zratio1[x, ], zratio2 = zratio2[x, ], tol = tol)})
  }
  return(out)
}

r_sol <- function(type1, type2, tau, zratio1, zratio2, tol) {
  f <- function(r)(bridge(type1 = type1, type2 = type2, r = r, zratio1 = zratio1, zratio2 = zratio2) - tau)^2
  op <- tryCatch(optimize(f, lower = -0.999, upper = 0.999, tol = tol)[1], error = function(e) 100)
  if(op == 100) {
    warning("Optimize returned error one of the pairwise correlations, returning NA")
    out <- NA
  } else {
    out <- unlist(op)
  }
}

bridge <- function(type1, type2, r, zratio1, zratio2) {
  if (type1 == "binary" & type2 == "continuous") {
    out <- bridgeF_bc(r = r, zratio1 = zratio1)
  } else if (type1 == "binary" & type2 == "binary") {
    out <- bridgeF_bb(r = r, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "trunc" & type2 == "continuous") {
    out <- bridgeF_tc(r = r, zratio1 = zratio1)
  } else if (type1 == "trunc" & type2 == "binary") {
    out <- bridgeF_tb(r = r, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "trunc" & type2 == "trunc") {
    out <- bridgeF_tt(r = r, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "ternary" & type2 == "continuous") {
    out <- bridgeF_nc(r = r, zratio1 = zratio1)
  } else if (type1 == "ternary" & type2 == "binary") {
    out <- bridgeF_nb(r = r, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "ternary" & type2 == "trunc") {
    out <- bridgeF_nt(r = r, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "ternary" & type2 == "ternary") {
    out <- bridgeF_nn(r = r, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "dtrunc" & type2 == "continuous") {
    out <- bridgeF_dc(r = r, zratio1 = zratio1)
  } else {
    stop("Unrecognized type of variables. Should be one of continuous, binary, trunc or ternary.")
  }
  return(out)
}
bridgeF_bc <- function(r, zratio1){
  # binary and continuous
  de1 <- stats::qnorm(zratio1)
  res <- as.numeric( 4*fMultivar::pnorm2d(de1, 0, rho = r/sqrt(2)) - 2*zratio1 )
  return(res)
}

bridgeF_bb <- function(r, zratio1, zratio2){
  # binary and binary
  de1 <- stats::qnorm(zratio1)
  de2 <- stats::qnorm(zratio2)
  res <- as.numeric(2 * (fMultivar::pnorm2d(de1, de2, rho = r) - zratio1*zratio2))
  return(res)
}

bridgeF_tc <- function(r, zratio1){
  # truncated and continuous
  de1 <- stats::qnorm(zratio1)
  mat2 <- matrix(c(1, 1/sqrt(2), r/sqrt(2),
                   1/sqrt(2), 1, r,
                   r/sqrt(2), r, 1), nrow = 3)
  res <- as.numeric(-2 * fMultivar::pnorm2d(-de1, 0, rho = 1/sqrt(2)) +
                      4 * mnormt::pmnorm(c(-de1, 0, 0), mean = rep(0, 3), varcov = mat2))
  return(res)
}

bridgeF_tb <- function(r, zratio1, zratio2){
  # truncated and binary
  de1 <- stats::qnorm(zratio1)
  de2 <- stats::qnorm(zratio2)
  mat1 <- matrix(c(1, -r, 1/sqrt(2),
                   -r, 1, -r/sqrt(2),
                   1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  mat2 <- matrix(c(1, 0, -1/sqrt(2),
                   0, 1, -r/sqrt(2),
                   -1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res <- as.numeric(2 * (1-zratio1) * (zratio2) -
      2 * mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat1) -
      2 * mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat2))
  return(res)
}

bridgeF_tt <- function(r, zratio1, zratio2){
  # truncated and truncated
  de1 <- stats::qnorm(zratio1)
  de2 <- stats::qnorm(zratio2)

  mat1 <- matrix(c(1, 0, 1/sqrt(2), -r/sqrt(2),
                   0, 1, -r/sqrt(2), 1/sqrt(2),
                   1/sqrt(2), -r/sqrt(2), 1, -r,
                   -r/sqrt(2), 1/sqrt(2), -r, 1), nrow = 4)
  mat2 <- matrix(c(1, r, 1/sqrt(2), r/sqrt(2),
                   r, 1, r/sqrt(2), 1/sqrt(2),
                   1/sqrt(2), r/sqrt(2), 1, r,
                   r/sqrt(2), 1/sqrt(2), r, 1), nrow = 4)

  res <- as.numeric(-2 * mnormt::pmnorm(c(-de1, -de2, 0, 0), mean = rep(0, 4), varcov = mat1) +
                    2 * mnormt::pmnorm(c(-de1, -de2, 0, 0), mean = rep(0, 4), varcov = mat2))
  return(res)
}

bridgeF_nc <- function(r, zratio1){
  # ternary and continuous
  de1 <- stats::qnorm(zratio1)
  mat <- matrix(c(1, 0, r/sqrt(2),
                   0, 1, -r/sqrt(2),
                   r/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res <- as.numeric(4 * fMultivar::pnorm2d(de1[2], 0, rho = r/sqrt(2)) - 2 * zratio1[2] +
         4 * mnormt::pmnorm(c(de1[1], de1[2], 0), mean = rep(0, 3), varcov = mat) -
         2 * zratio1[1]*zratio1[2])
  return(res)
}

bridgeF_nb <- function(r, zratio1, zratio2){
  # ternary and binary
  de1 <- stats::qnorm(zratio1)
  de2 <- stats::qnorm(zratio2)

  res <- as.numeric(2 * fMultivar::pnorm2d(de2, de1[2], rho = r) * (1 - zratio1[1]) -
                      2 * zratio1[2] * (zratio2 - fMultivar::pnorm2d(de2, de1[1], rho = r)))
  return(res)
}

bridgeF_nt <- function(r, zratio1, zratio2){
  # ternary and binary
  de1 <- stats::qnorm(zratio1)
  de2 <- stats::qnorm(zratio2)

  mat1 <- matrix(c(1, 0, 0,
                   0, 1, r,
                   0, r, 1), nrow = 3)

  mat2 <- matrix(c(1, 0, 0, r/sqrt(2),
                   0, 1, -r, r/sqrt(2),
                   0, -r, 1, -1/sqrt(2),
                   r/sqrt(2), r/sqrt(2), -1/sqrt(2), 1), nrow = 4)

  mat3 <- matrix(c(1, 0, r, r/sqrt(2),
                   0, 1, 0, r/sqrt(2),
                   r, 0, 1, 1/sqrt(2),
                   r/sqrt(2), r/sqrt(2), 1/sqrt(2), 1), nrow = 4)

  res <- as.numeric(- 2 * (1 - zratio1[1]) * zratio1[2]
                    + 2 * mnormt::pmnorm(c(-de1[1], de1[2], de2), mean = rep(0, 3), varcov = mat1)
                    + 2 * mnormt::pmnorm(c(-de1[1], de1[2], -de2, 0), mean = rep(0, 4), varcov = mat2)
                    + 2 * mnormt::pmnorm(c(-de1[1], de1[2], -de2, 0), mean = rep(0, 4), varcov = mat3))
  return(res)
}

bridgeF_nn <- function(r, zratio1, zratio2){
  # ternary and ternary
  de1 <- stats::qnorm(zratio1)
  de2 <- stats::qnorm(zratio2)

  res <- as.numeric(2 * fMultivar::pnorm2d(de1[2], de2[2], rho = r) *
                      fMultivar::pnorm2d(-de1[1], -de2[1], rho = r) -
         2 * (zratio1[2] - fMultivar::pnorm2d(de1[2], de2[1], rho = r)) *
           (zratio2[2] - fMultivar::pnorm2d(de1[1], de2[2], rho = r)))
  return(res)
}

# bridgeF_dc <- function(r, zratio1){
#   de1 <- stats::qnorm(zratio1)
#
#   mat1 <- matrix(c(1, 0, r/sqrt(2),
#                    0, 1, r/sqrt(2),
#                    r/sqrt(2), r/sqrt(2), 1), nrow = 3)
#   mat2 <- matrix(c(1, 0, - 1/sqrt(2),
#                    0, 1, - 1/sqrt(2),
#                    - 1/sqrt(2), - 1/sqrt(2), 1), nrow = 3)
#   mat3 <- matrix(c(1, 0, - 1/sqrt(2), - r/sqrt(2),
#                    0, 1, - 1/sqrt(2), - r/sqrt(2),
#                    - 1/sqrt(2), - 1/sqrt(2), 1, r,
#                    - r/sqrt(2), - r/sqrt(2), r, 1), nrow = 4)
#
#   res <- as.numeric(2 * zratio1[1]^2 - 4 * zratio1[1] + 2 * zratio1[1] * zratio1[2] + 2 * zratio1[2]^2 - 2 * zratio1[2]
#                     + 4 * mnormt::pmnorm(c(de1[1], -de1[1], 0), mean = rep(0, 3), varcov = mat1)
#                     + 4 * mnormt::pmnorm(c(de1[2], -de1[2], 0), mean = rep(0, 3), varcov = mat1)
#                     + 4 * mnormt::pmnorm(c(de1[1], -de1[2], 0), mean = rep(0, 3), varcov = mat1)
#                     - 2 * mnormt::pmnorm(c(-de1[1], de1[2], 0), mean = rep(0, 3), varcov = mat2)
#                     + 4 * mnormt::pmnorm(c(-de1[1], de1[2], 0, 0), mean = rep(0, 4), varcov = mat3))
#   return(res)
# }
############################################################################################
# For multilinear interpolation approximation for bridge Inverse
############################################################################################
bound_bc <- function(zratio1){2 * zratio1 * (1 - zratio1)}
bound_bb <- function(zratio1, zratio2){2 * pmin(zratio1, zratio2)*(1-pmax(zratio1, zratio2))}
bound_tc <- function(zratio1){1 - zratio1^2}
bound_tb <- function(zratio1, zratio2){2 * pmax(zratio2, 1 - zratio2) * (1 - pmax(zratio2, 1 - zratio2, zratio1))}
bound_tt <- function(zratio1, zratio2){1 - pmax(zratio1, zratio2)^2}
bound_nc <- function(zratio1){2 * (zratio1[ , 1] * (1 - zratio1[ , 1]) + (1 - zratio1[ , 2]) * (zratio1[ , 2] - zratio1[ , 1]))}
bound_nb <- function(zratio1, zratio2){2 * pmin(zratio1[ , 1] * (1 - zratio1[ , 1]) + (1 - zratio1[ , 2]) * (zratio1[ , 2] - zratio1[ , 1]),
                                                zratio2 * (1 - zratio2))}
bound_nt <- function(zratio1, zratio2){1 - pmax(zratio1[ , 1], zratio1[ , 2] - zratio1[ , 1], 1 - zratio1[ , 2], zratio2)^2}
bound_nn <- function(zratio1, zratio2){2 * pmin(zratio1[ , 1] * (1 - zratio1[ , 1]) + (1 - zratio1[ , 2]) * (zratio1[ , 2] - zratio1[ , 1]),
                                                zratio2[ , 1] * (1 - zratio2[ , 1]) + (1 - zratio2[ , 2]) * (zratio2[ , 2] - zratio2[ , 1]))}
############################################################################################

cutoff <- function(type1, type2, tau, zratio1, zratio2, method, ratio){
  if (method == "original") {
    out <- rep(TRUE, length(tau))
  } else if (method == "ml") {
    out <- rep(FALSE, length(tau))
  } else if (type1 == "binary" & type2 == "continuous") {
    out <- tau > ratio * bound_bc(zratio1 = zratio1)
  } else if (type1 == "binary" & type2 == "binary") {
    out <- tau > ratio * bound_bb(zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "trunc" & type2 == "continuous") {
    out <- tau > ratio * bound_tc(zratio1 = zratio1)
  } else if (type1 == "trunc" & type2 == "binary") {
    out <- tau > ratio * bound_tb(zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "trunc" & type2 == "trunc") {
    out <- tau > ratio * bound_tt(zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "ternary" & type2 == "continuous") {
    out <- tau > ratio * bound_nc(zratio1 = zratio1)
  } else if (type1 == "ternary" & type2 == "binary") {
    out <- tau > ratio * bound_nb(zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "ternary" & type2 == "trunc") {
    out <- tau > ratio * bound_nt(zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "ternary" & type2 == "ternary") {
    out <- tau > ratio * bound_nn(zratio1 = zratio1, zratio2 = zratio2)
  # } else if (type1 == "dtrunc" & type2 == "continuous") {
  #   out <- rep(TRUE, length(tau))
  } else {
    stop("Unrecognized type of variables. Should be one of continuous, binary or trunc.")
  }
  return(out)
}



############################################################################################
# Select which bridge inverse function based on the combination of variable types
############################################################################################

r_ml <- function(type1, type2, tau, zratio1, zratio2) {
  if (type1 == "binary" & type2 == "continuous") {
    out <- BCipol(t(cbind(tau / bound_bc(zratio1 = zratio1), zratio1))) / 10^7
  } else if (type1 == "binary" & type2 == "binary") {
    out <- BBipol(t(cbind(tau / bound_bb(zratio1 = zratio1, zratio2 = zratio2), zratio1, zratio2))) / 10^7
  } else if (type1 == "trunc" & type2 == "continuous") {
    out <- TCipol(t(cbind(tau / bound_tc(zratio1 = zratio1), zratio1))) / 10^7
  } else if (type1 == "trunc" & type2 == "binary") {
    out <- TBipol(t(cbind(tau / bound_tb(zratio1 = zratio1, zratio2 = zratio2), zratio1, zratio2))) / 10^7
  } else if (type1 == "trunc" & type2 == "trunc") {
    out <- TTipol(t(cbind(tau / bound_tt(zratio1 = zratio1, zratio2 = zratio2), zratio1, zratio2))) / 10^7
  } else if (type1 == "ternary" & type2 == "continuous") {
    out <- NCipol(t(cbind(tau / bound_nc(zratio1 = zratio1), zratio1[ , 1] / zratio1[ , 2], zratio1[ , 2]))) / 10^7
  } else if (type1 == "ternary" & type2 == "binary") {
    out <- NBipol(t(cbind(tau / bound_nb(zratio1 = zratio1, zratio2 = zratio2), zratio1[ , 1] / zratio1[ , 2], zratio1[ , 2], zratio2))) / 10^7
  } else if (type1 == "ternary" & type2 == "trunc") {
    out <- NTipol(t(cbind(tau / bound_nt(zratio1 = zratio1, zratio2 = zratio2), zratio1[ , 1] / zratio1[ , 2], zratio1[ , 2], zratio2))) / 10^7
  } else if (type1 == "ternary" & type2 == "ternary") {
    out <- NNipol(t(cbind(tau / bound_nn(zratio1 = zratio1, zratio2 = zratio2), zratio1[ , 1] / zratio1[ , 2], zratio1[ , 2], zratio2[ , 1] / zratio2[ , 2], zratio2[ , 2]))) / 10^7
  } else {
    stop("Unrecognized type of variables. Should be one of continuous, binary or trunc.")
  }
  return(out)
}
