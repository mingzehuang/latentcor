#' @importFrom stats qnorm
#' @importFrom mnormt pmnorm
#' @importFrom fMultivar pnorm2d
#'
NULL

copula_list = vector(mode = "list", length = 2)
names(copula_list) = c("exp", "cube")
copula_list[[1]] = function(z) {exp(z)}
copula_list[[2]] = function(z) {z^3}

type_list = vector(mode = "list", length = 4)
names(type_list) = c("con", "bin", "tru", "ter")
for (i in seq(length(type_list))) {type_list[[i]] = i - 1}

transform_list = zratio_list = vector(mode = "list", length = 3)

transform_list[[1]] = function(u, pi) {
  if(length(pi) != 1) {
    stop("Proportion of zeros should be specified by a scalar pi for binary data.")
  } else {
    q = quantile(u, pi); x = ifelse(u > q, 1, 0)
  }
  return(x)
}

transform_list[[2]] = function(u, pi) {
  if(length(pi) != 1) {
    stop("Proportion of zeros should be specified by a scalar pi for truncated data.")
  } else {
    q = quantile(u, pi); x = ifelse(u > q, u, q) - q
  }
  return(x)
}

transform_list[[3]] = function(u, pi) {
  if(length(pi) != 2) {
    stop("Proportion of zeros and ones should be specified by a vector pi for ternary data.")
  } else {
  q = quantile(u, pi); x = rep(1, length(u)); x[u > q[2]] = 2; x[u <= q[1]] = 0
  }
  return(x)
}

zratio_list[[1]] = function(X) {colMeans(as.matrix(X) == 0)}
zratio_list[[2]] = function(X) {colMeans(as.matrix(X) == 0)}
zratio_list[[3]] = function(X) {mattolist(rbind(colMeans(as.matrix(X) == 0), 1 - colMeans(as.matrix(X) == 2)))}

bridge_list = bound_list = vector(mode = "list", length = 9)
names(bridge_list) = names(bound_list) = c("10", "11", "20", "21", "22", "30", "31", "32", "33")

bridge_list[[1]] = function(r, zratio1){
  # binary and continuous
  de1 = stats::qnorm(zratio1)
  res = as.numeric(4*fMultivar::pnorm2d(de1, 0, rho = r/sqrt(2)) - 2*zratio1)
  return(res)
}

bridge_list[[2]] = function(r, zratio1, zratio2){
  # binary and binary
  de1 = stats::qnorm(zratio1); de2 = stats::qnorm(zratio2)
  res = as.numeric(2 * (fMultivar::pnorm2d(de1, de2, rho = r) - zratio1*zratio2))
  return(res)
}

bridge_list[[3]] = function(r, zratio1){
  # truncated and continuous
  de1 = stats::qnorm(zratio1)
  mat2 = matrix(c(1, 1/sqrt(2), r/sqrt(2), 1/sqrt(2), 1, r, r/sqrt(2), r, 1), nrow = 3)
  res = as.numeric(-2 * fMultivar::pnorm2d(-de1, 0, rho = 1/sqrt(2)) + 4 * mnormt::pmnorm(c(-de1, 0, 0), mean = rep(0, 3), varcov = mat2))
  return(res)
}

bridge_list[[4]] = function(r, zratio1, zratio2){
  # truncated and binary
  de1 = stats::qnorm(zratio1); de2 = stats::qnorm(zratio2)
  mat1 = matrix(c(1, -r, 1/sqrt(2), -r, 1, -r/sqrt(2), 1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  mat2 = matrix(c(1, 0, -1/sqrt(2), 0, 1, -r/sqrt(2), -1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res = as.numeric(2 * (1-zratio1) * (zratio2) - 2 * mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat1)
                   - 2 * mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat2))
  return(res)
}

bridge_list[[5]] = function(r, zratio1, zratio2){
  # truncated and truncated
  de1 = stats::qnorm(zratio1); de2 = stats::qnorm(zratio2)
  mat1 = matrix(c(1, 0, 1/sqrt(2), -r/sqrt(2), 0, 1, -r/sqrt(2), 1/sqrt(2), 1/sqrt(2), -r/sqrt(2), 1, -r, -r/sqrt(2), 1/sqrt(2), -r, 1), nrow = 4)
  mat2 = matrix(c(1, r, 1/sqrt(2), r/sqrt(2), r, 1, r/sqrt(2), 1/sqrt(2), 1/sqrt(2), r/sqrt(2), 1, r, r/sqrt(2), 1/sqrt(2), r, 1), nrow = 4)
  res = as.numeric(-2 * mnormt::pmnorm(c(-de1, -de2, 0, 0), mean = rep(0, 4), varcov = mat1)
                   + 2 * mnormt::pmnorm(c(-de1, -de2, 0, 0), mean = rep(0, 4), varcov = mat2))
  return(res)
}

bridge_list[[6]] = function(r, zratio1){
  # ternary and continuous
  de1 = stats::qnorm(zratio1)
  mat = matrix(c(1, 0, r/sqrt(2), 0, 1, -r/sqrt(2), r/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res = as.numeric(4 * fMultivar::pnorm2d(de1[2], 0, rho = r/sqrt(2)) - 2 * zratio1[2]
                   + 4 * mnormt::pmnorm(c(de1[1], de1[2], 0), mean = rep(0, 3), varcov = mat) - 2 * zratio1[1]*zratio1[2])
  return(res)
}

bridge_list[[7]] = function(r, zratio1, zratio2){
  # ternary and binary
  de1 = stats::qnorm(zratio1); de2 = stats::qnorm(zratio2)
  res = as.numeric(2 * fMultivar::pnorm2d(de2, de1[2], rho = r) * (1 - zratio1[1])
                   - 2 * zratio1[2] * (zratio2 - fMultivar::pnorm2d(de2, de1[1], rho = r)))
  return(res)
}

bridge_list[[8]] = function(r, zratio1, zratio2){
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

bridge_list[[9]] = function(r, zratio1, zratio2){
  # ternary and ternary
  de1 = stats::qnorm(zratio1); de2 = stats::qnorm(zratio2)
  res = as.numeric(2 * fMultivar::pnorm2d(de1[2], de2[2], rho = r) * fMultivar::pnorm2d(-de1[1], -de2[1], rho = r)
                   - 2 * (zratio1[2] - fMultivar::pnorm2d(de1[2], de2[1], rho = r)) * (zratio2[2] - fMultivar::pnorm2d(de1[1], de2[2], rho = r)))
  return(res)
}

bound_list[[1]] = function(zratio1){2 * zratio1 * (1 - zratio1)}
bound_list[[2]] = function(zratio1, zratio2){2 * pmin(zratio1, zratio2)*(1-pmax(zratio1, zratio2))}
bound_list[[3]] = function(zratio1){1 - zratio1^2}
bound_list[[4]] = function(zratio1, zratio2){2 * pmax(zratio2, 1 - zratio2) * (1 - pmax(zratio2, 1 - zratio2, zratio1))}
bound_list[[5]] = function(zratio1, zratio2){1 - pmax(zratio1, zratio2)^2}
bound_list[[6]] = function(zratio1){2 * (zratio1[1] * (1 - zratio1[1]) + (1 - zratio1[2]) * (zratio1[2] - zratio1[1]))}
bound_list[[7]] = function(zratio1, zratio2){2 * pmin(zratio1[1] * (1 - zratio1[1]) + (1 - zratio1[2]) * (zratio1[2] - zratio1[1]), zratio2 * (1 - zratio2))}
bound_list[[8]] = function(zratio1, zratio2){1 - pmax(zratio1[1], zratio1[2] - zratio1[1], 1 - zratio1[2], zratio2)^2}
bound_list[[9]] = function(zratio1, zratio2){2 * pmin(zratio1[1] * (1 - zratio1[1]) + (1 - zratio1[2]) * (zratio1[2] - zratio1[1]),
                                                zratio2[1] * (1 - zratio2[1]) + (1 - zratio2[2]) * (zratio2[2] - zratio2[1]))}
