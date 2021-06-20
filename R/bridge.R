#' @importFrom stats qnorm
#' @importFrom mnormt pmnorm
#' @importFrom fMultivar pnorm2d
#'
NULL

bridge_10 = function(r, zratio1){
  # binary and continuous
  de1 = stats::qnorm(zratio1)
  res = as.numeric(4*fMultivar::pnorm2d(de1, 0, rho = r/sqrt(2)) - 2*zratio1)
  return(res)
}

bridge_11 = function(r, zratio1, zratio2){
  # binary and binary
  de1 = stats::qnorm(zratio1); de2 = stats::qnorm(zratio2)
  res = as.numeric(2 * (fMultivar::pnorm2d(de1, de2, rho = r) - zratio1*zratio2))
  return(res)
}

bridge_20 = function(r, zratio1){
  # truncated and continuous
  de1 = stats::qnorm(zratio1)
  mat2 = matrix(c(1, 1/sqrt(2), r/sqrt(2), 1/sqrt(2), 1, r, r/sqrt(2), r, 1), nrow = 3)
  res = as.numeric(-2 * fMultivar::pnorm2d(-de1, 0, rho = 1/sqrt(2)) + 4 * mnormt::pmnorm(c(-de1, 0, 0), mean = rep(0, 3), varcov = mat2))
  return(res)
}

bridge_21 = function(r, zratio1, zratio2){
  # truncated and binary
  de1 = stats::qnorm(zratio1); de2 = stats::qnorm(zratio2)
  mat1 = matrix(c(1, -r, 1/sqrt(2), -r, 1, -r/sqrt(2), 1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  mat2 = matrix(c(1, 0, -1/sqrt(2), 0, 1, -r/sqrt(2), -1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res = as.numeric(2 * (1-zratio1) * (zratio2) - 2 * mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat1)
                   - 2 * mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat2))
  return(res)
}

bridge_22 = function(r, zratio1, zratio2){
  # truncated and truncated
  de1 = stats::qnorm(zratio1); de2 = stats::qnorm(zratio2)
  mat1 = matrix(c(1, 0, 1/sqrt(2), -r/sqrt(2), 0, 1, -r/sqrt(2), 1/sqrt(2), 1/sqrt(2), -r/sqrt(2), 1, -r, -r/sqrt(2), 1/sqrt(2), -r, 1), nrow = 4)
  mat2 = matrix(c(1, r, 1/sqrt(2), r/sqrt(2), r, 1, r/sqrt(2), 1/sqrt(2), 1/sqrt(2), r/sqrt(2), 1, r, r/sqrt(2), 1/sqrt(2), r, 1), nrow = 4)
  res = as.numeric(-2 * mnormt::pmnorm(c(-de1, -de2, 0, 0), mean = rep(0, 4), varcov = mat1)
                   + 2 * mnormt::pmnorm(c(-de1, -de2, 0, 0), mean = rep(0, 4), varcov = mat2))
  return(res)
}

bridge_30 = function(r, zratio1){
  # ternary and continuous
  de1 = stats::qnorm(zratio1)
  mat = matrix(c(1, 0, r/sqrt(2), 0, 1, -r/sqrt(2), r/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res = as.numeric(4 * fMultivar::pnorm2d(de1[2], 0, rho = r/sqrt(2)) - 2 * zratio1[2]
                   + 4 * mnormt::pmnorm(c(de1[1], de1[2], 0), mean = rep(0, 3), varcov = mat) - 2 * zratio1[1]*zratio1[2])
  return(res)
}

bridge_31 = function(r, zratio1, zratio2){
  # ternary and binary
  de1 = stats::qnorm(zratio1); de2 = stats::qnorm(zratio2)
  res = as.numeric(2 * fMultivar::pnorm2d(de2, de1[2], rho = r) * (1 - zratio1[1])
                   - 2 * zratio1[2] * (zratio2 - fMultivar::pnorm2d(de2, de1[1], rho = r)))
  return(res)
}

bridge_32 = function(r, zratio1, zratio2){
  # ternary and binary
  de1 = stats::qnorm(zratio1); de2 = stats::qnorm(zratio2)
  mat1 = matrix(c(1, 0, 0, 0, 1, r, 0, r, 1), nrow = 3)
  mat2 = matrix(c(1, 0, 0, r/sqrt(2), 0, 1, -r, r/sqrt(2), 0, -r, 1, -1/sqrt(2), r/sqrt(2), r/sqrt(2), -1/sqrt(2), 1), nrow = 4)
  mat3 = matrix(c(1, 0, r, r/sqrt(2), 0, 1, 0, r/sqrt(2), r, 0, 1, 1/sqrt(2), r/sqrt(2), r/sqrt(2), 1/sqrt(2), 1), nrow = 4)
  res = as.numeric(- 2 * (1 - zratio1[1]) * zratio1[2]
                    + 2 * mnormt::pmnorm(c(-de1[1], de1[2], de2), mean = rep(0, 3), varcov = mat1)
                    + 2 * mnormt::pmnorm(c(-de1[1], de1[2], -de2, 0), mean = rep(0, 4), varcov = mat2)
                    + 2 * mnormt::pmnorm(c(-de1[1], de1[2], -de2, 0), mean = rep(0, 4), varcov = mat3))
  return(res)
}

bridge_33 = function(r, zratio1, zratio2){
  # ternary and ternary
  de1 = stats::qnorm(zratio1); de2 = stats::qnorm(zratio2)
  res = as.numeric(2 * fMultivar::pnorm2d(de1[2], de2[2], rho = r) * fMultivar::pnorm2d(-de1[1], -de2[1], rho = r)
                   - 2 * (zratio1[2] - fMultivar::pnorm2d(de1[2], de2[1], rho = r)) * (zratio2[2] - fMultivar::pnorm2d(de1[1], de2[2], rho = r)))
  return(res)
}

bound_10 = function(zratio1){2 * zratio1 * (1 - zratio1)}
bound_11 = function(zratio1, zratio2){2 * pmin(zratio1, zratio2)*(1-pmax(zratio1, zratio2))}
bound_20 = function(zratio1){1 - zratio1^2}
bound_21 = function(zratio1, zratio2){2 * pmax(zratio2, 1 - zratio2) * (1 - pmax(zratio2, 1 - zratio2, zratio1))}
bound_22 = function(zratio1, zratio2){1 - pmax(zratio1, zratio2)^2}
bound_30 = function(zratio1){2 * (zratio1[1] * (1 - zratio1[1]) + (1 - zratio1[2]) * (zratio1[2] - zratio1[1]))}
bound_31 = function(zratio1, zratio2){2 * pmin(zratio1[1] * (1 - zratio1[1]) + (1 - zratio1[2]) * (zratio1[2] - zratio1[1]), zratio2 * (1 - zratio2))}
bound_32 = function(zratio1, zratio2){1 - pmax(zratio1[1], zratio1[2] - zratio1[1], 1 - zratio1[2], zratio2)^2}
bound_33 = function(zratio1, zratio2){2 * pmin(zratio1[1] * (1 - zratio1[1]) + (1 - zratio1[2]) * (zratio1[2] - zratio1[1]),
                                                zratio2[1] * (1 - zratio2[1]) + (1 - zratio2[2]) * (zratio2[2] - zratio2[1]))}
