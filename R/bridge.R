#'
#' @importFrom chebpol ipol
#'
#'
NULL

r_sol <- function(type1, type2, tau, zratio1, zratio2, method, tol = NULL) {
  if (method == "original") {
    f <- function(r)(bridge(type1 = type1, type2 = type2, r = r, zratio1 = zratio1, zratio2 = zratio2) - tau)^2
    op <- tryCatch(optimize(f, lower = -0.999, upper = 0.999, tol = tol)[1], error = function(e) 100)
      if(op == 100) {
        warning("Optimize returned error one of the pairwise correlations, returning NA")
        out <- NA
    } else {
      out <- unlist(op)
    }
  } else if (method == "ml") {
    out <- bridgeInv(type1 = type1, type2 = type2, tau = tau, zratio1 = zratio1, zratio2 = zratio2)
  }
  return(out)
}

bridge <- function(type1, type2, r, zratio1, zratio2) {
  if (type1 == "binary" & type2 == "continuous") {
    out <- bridgeF_bc(r = r, zratio1 = zratio1)
  } else if (type1 == "continuous" & type2 == "binary") {
    out <- bridgeF_bc(r = r, zratio1 = zratio2)
  } else if (type1 == "binary" & type2 == "binary") {
    out <- bridgeF_bb(r = r, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "trunc" & type2 == "continuous") {
    out <- bridgeF_tc(r = r, zratio1 = zratio1)
  } else if (type1 == "continuous" & type2 == "trunc") {
    out <- bridgeF_tc(r = r, zratio1 = zratio2)
  } else if (type1 == "trunc" & type2 == "binary") {
    out <- bridgeF_tb(r = r, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "binary" & type2 == "trunc") {
    out <- bridgeF_tb(r = r, zratio1 = zratio2, zratio2 = zratio1)
  } else if (type1 == "trunc" & type2 == "trunc") {
    out <- bridgeF_tt(r = r, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "ternary" & type2 == "continuous") {
    out <- bridgeF_nc(r = r, zratio1 = zratio1)
  } else if (type1 == "continuous" & type2 == "ternary") {
    out <- bridgeF_nc(r = r, zratio1 = zratio2)
  } else if (type1 == "ternary" & type2 == "binary") {
    out <- bridgeF_nb(r = r, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "binary" & type2 == "ternary") {
    out <- bridgeF_nb(r = r, zratio1 = zratio2, zratio2 = zratio1)
  } else if (type1 == "ternary" & type2 == "ternary") {
    out <- bridgeF_nn(r = r, zratio1 = zratio1, zratio2 = zratio2)
  # } else if (type1 == "ternary" & type2 == "trunc") {bridge_select <- bridgeF_nt
  # } else if (type1 == "trunc" & type2 == "ternary") {bridge_select <- bridgeF_tn
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
  de1 <- stats::qnorm(zratio1[1])
  de2 <- stats::qnorm(zratio1[2])
  mat <- matrix(c(1, 0, r/sqrt(2),
                   0, 1, -r/sqrt(2),
                   r/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res <- as.numeric(4 * fMultivar::pnorm2d(de2, 0, rho = r/sqrt(2)) - 2 * zratio1[2] +
         4 * mnormt::pmnorm(c(de1, de2, 0), mean = rep(0, 3), varcov = mat) -
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

############################################################################################
# For multilinear interpolation approximation for bridge Inverse
############################################################################################

############################################################################################
# Cutoff criteria based on the combination of variable types
############################################################################################
cutoff_bc <- function(zratio1){0.9 * 2 * zratio1 * (1 - zratio1)}
cutoff_bb <- function(zratio1, zratio2){0.9 * 2 * pmin(zratio1, zratio2)*(1-pmax(zratio1, zratio2))}
cutoff_tc <- function(zratio1){0.9 * (1 - zratio1^2)}
cutoff_tb <- function(zratio1, zratio2){0.9 * 2 * pmax(zratio2, 1 - zratio2) * (1 - pmax(zratio2, 1 - zratio2, zratio1))}
cutoff_tt <- function(zratio1, zratio2){0.9 * (1 - pmax(zratio1, zratio2)^2)}
cutoff_nc <- function(zratio1){0.9 * 2 * (zratio1[ , 1] * (zratio1[ , 2] - zratio1[ , 1]) + (1 - zratio1[ , 2]) * zratio1[ , 2])}
cutoff_nb <- function(zratio1, zratio2){0.9 * 2 * pmin(zratio1[ , 1] * (zratio1[ , 2] - zratio1[ , 1]) + (1 - zratio1[ , 2]) * zratio1[ , 2], zratio2 * (1 - zratio2))}
cutoff_nn <- function(zratio1, zratio2){0.9 * 2 * pmin(zratio1[ , 1] * (zratio1[ , 2] - zratio1[ , 1]) + (1 - zratio1[ , 2]) * zratio1[ , 2],
                                                       zratio2[ , 1] * (zratio2[ , 2] - zratio2[ , 1]) + (1 - zratio2[ , 2]) * zratio2[ , 2])}

cutoff <- function(type1, type2, tau, zratio1, zratio2, method){
  if (method == "original") {
    out <- rep(TRUE, length(tau))
  } else if (method == "ml") {
    out <- rep(FALSE, length(tau))
  } else if (type1 == "binary" & type2 == "continuous") {
    out <- c(tau) > cutoff_bc(zratio1 = zratio1)
  } else if (type1 == "continuous" & type2 == "binary") {
    out <- c(tau) > cutoff_bc(zratio1 = zratio2)
  } else if (type1 == "binary" & type2 == "binary") {
    out <- c(tau) > cutoff_bb(zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "trunc" & type2 == "continuous") {
    out <- c(tau) > cutoff_tc(zratio1 = zratio1)
  } else if (type1 == "continuous" & type2 == "trunc") {
    out <- c(tau) > cutoff_tc(zratio1 = zratio2)
  } else if (type1 == "trunc" & type2 == "binary") {
    out <- c(tau) > cutoff_tb(zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "binary" & type2 == "trunc") {
    out <- c(tau) > cutoff_tb(zratio1 = zratio2, zratio2 = zratio1)
  } else if (type1 == "trunc" & type2 == "trunc") {
    out <- c(tau) > cutoff_tt(zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "ternary" & type2 == "continuous") {
    out <- c(tau) > cutoff_nc(zratio1 = zratio1)
  } else if (type1 == "continuous" & type2 == "ternary") {
    out <- c(tau) > cutoff_nc(zratio1 = zratio2)
  } else if (type1 == "ternary" & type2 == "binary") {
    out <- c(tau) > cutoff_nb(zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "binary" & type2 == "ternary") {
    out <- c(tau) > cutoff_nb(zratio1 = zratio2, zratio2 = zratio1)
  } else if (type1 == "ternary" & type2 == "ternary") {
    out <- c(tau) > cutoff_nn(zratio1 = zratio1, zratio2 = zratio2)
  } else {
    stop("Unrecognized type of variables. Should be one of continuous, binary or trunc.")
  }
  return(out)
}



############################################################################################
# Select which bridge inverse function based on the combination of variable types
############################################################################################

bridgeInv <- function(type1, type2, tau, zratio1, zratio2) {
  if (type1 == "binary" & type2 == "continuous") {
    out <- bridgeInv_bc(tau = tau, zratio1 = zratio1)
  } else if (type1 == "continuous" & type2 == "binary") {
    out <- bridgeInv_bc(tau = tau, zratio1 = zratio2)
  } else if (type1 == "binary" & type2 == "binary") {
    out <- bridgeInv_bb(tau = tau, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "trunc" & type2 == "continuous") {
    out <- bridgeInv_tc(tau = tau, zratio1 = zratio1)
  } else if (type1 == "continuous" & type2 == "trunc") {
    out <- bridgeInv_tc(tau = tau, zratio1 = zratio2)
  } else if (type1 == "trunc" & type2 == "trunc") {
    out <- bridgeInv_tt(tau = tau, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "trunc" & type2 == "binary") {
    out <- bridgeInv_tb(tau = tau, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "binary" & type2 == "trunc") {
    out <- bridgeInv_tb(tau = tau, zratio1 = zratio2, zratio2 = zratio1)
  } else if (type1 == "ternary" & type2 == "binary") {
    out <- bridgeInv_nb(tau = tau, zratio1 = zratio1, zratio2 = zratio2)
  } else if (type1 == "binary" & type2 == "ternary") {
    out <- bridgeInv_nb(tau = tau, zratio1 = zratio2, zratio2 = zratio1)
  } else if (type1 == "ternary" & type2 == "continuous") {
    out <- bridgeInv_nc(tau = tau, zratio1 = zratio1)
  } else if (type1 == "continuous" & type2 == "ternary") {
    out <- bridgeInv_nc(tau = tau, zratio1 = zratio2)
  } else if (type1 == "ternary" & type2 == "ternary") {
    out <- bridgeInv_nn(tau = tau, zratio1 = zratio1, zratio2 = zratio2)
  } else {
    stop("Unrecognized type of variables. Should be one of continuous, binary or trunc.")
  }
  return(out)
}


# wrapper function for BC
bridgeInv_bc <- function(tau, zratio1){
  out <- BCipol(rbind(t(tau / (2 * zratio1 * (1 - zratio1))), t(zratio1))) / 10^7
  return(out)
}

# wrapper function
bridgeInv_bb <- function(tau, zratio1, zratio2){
  out <- BBipol(rbind(t(tau / (2 * pmin(zratio1, zratio2)*(1-pmax(zratio1, zratio2)))), t(zratio1), t(zratio2))) / 10^7
  return(out)
}

# wrapper functions
bridgeInv_tc <- function(tau, zratio1){
  out <- TCipol(rbind(t(tau / (1 - zratio1^2)), t(zratio1))) / 10^7
  return(out)
}

# wrapper functions
bridgeInv_tb <- function(tau, zratio1, zratio2){
  out <- TBipol(rbind(t(tau / (2 * pmax(zratio2, 1 - zratio2) * (1 - pmax(zratio2, 1 - zratio2, zratio1)))), t(zratio1), t(zratio2))) / 10^7
  return(out)
}

# wrapper function
bridgeInv_tt <- function(tau, zratio1, zratio2){
  out <- TTipol(rbind(t(tau / (1 - pmax(zratio1, zratio2)^2)), t(zratio1), t(zratio2))) / 10^7
  return(out)
}

# wrapper function
bridgeInv_nc <- function(tau, zratio1){
  out <- NCipol(rbind(t(tau / (2 * (zratio1[ , 1] * (zratio1[ , 2] - zratio1[ , 1]) + (1 - zratio1[ , 2]) * zratio1[ , 2]))), t(zratio1[ , 1] / zratio1[ , 2]), t(zratio1[ , 2]))) / 10^7
  return(out)
}

# wrapper function
bridgeInv_nb <- function(tau, zratio1, zratio2){
  out <- NBipol(rbind(t(tau / (2 * pmin(zratio1[ , 1] * (zratio1[ , 2] - zratio1[ , 1]) + (1 - zratio1[ , 2]) * zratio1[ , 2], zratio2 * (1 - zratio2)))), t(zratio1[ , 1] / zratio1[ , 2]), t(zratio1[ , 2]), t(zratio2))) / 10^7
  return(out)
}

# wrapper function
bridgeInv_nn <- function(tau, zratio1, zratio2){
  out <- NNipol(rbind(t(tau / (2 * pmin(zratio1[ , 1] * (zratio1[ , 2] - zratio1[ , 1]) + (1 - zratio1[ , 2]) * zratio1[ , 2],
                                        zratio2[ , 1] * (zratio2[ , 2] - zratio2[ , 1]) + (1 - zratio2[ , 2]) * zratio2[ , 2]))), t(zratio1[ , 1] / zratio1[ , 2]), t(zratio1[ , 2]), t(zratio2[ , 1] / zratio2[ , 2]), t(zratio2[ , 2]))) / 10^7
  return(out)
}

