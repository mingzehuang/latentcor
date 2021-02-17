bridge_select <- function(type1 = "trunc", type2 = "continuous") {
  if (type1 == "binary" & type2 == "binary") { bridge_select <- bridgeF_bb
  } else if (type1 == "trunc" & type2 == "trunc") { bridge_select <- bridgeF_tt
  } else if (type1 == "trunc" & type2 == "continuous") { bridge_select <- bridgeF_tc
  } else if (type1 == "continuous" & type2 == "trunc") { bridge_select <- bridgeF_ct
  } else if (type1 == "binary" & type2 == "continuous") { bridge_select <- bridgeF_bc
  } else if (type1 == "continuous" & type2 == "binary") { bridge_select <- bridgeF_cb
  } else if (type1 == "trunc" & type2 == "binary") { bridge_select <- bridgeF_tb
  } else if (type1 == "binary" & type2 == "trunc") { bridge_select <- bridgeF_bt
  } else if (type1 == "ternary" & type2 == "continuous") { bridge_select <- bridgeF_nc
  } else if (type1 == "ternary" & type2 == "ternary") { bridge_select <- bridgeF_nn
  } else if (type1 == "ordinal" & type2 == "continuous") { bridge_select <- bridgeF_oc
  } else {
    stop("Unrecognized type of variables. Should be one of continuous, binary or trunc.")
  }
}
bridgeF_bc <- function(r, zratio1, zratio2 = NULL){
  # binary and continuous
  de1 <- stats::qnorm(zratio1)
  res <- as.numeric( 4*fMultivar::pnorm2d(de1, 0, rho = r/sqrt(2)) - 2*zratio1 )
  return(res)
}
bridgeF_cb <- function(r, zratio1 = NULL, zratio2){
  # continuous and binary
  de2 <- stats::qnorm(zratio2)
  res <- as.numeric(4 * fMultivar::pnorm2d(0, de2, rho = r/sqrt(2)) - 2 * zratio2 )
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
bridgeF_bt <- function(r, zratio1, zratio2){
  # binary and truncated
  de1 <- stats::qnorm(zratio2)
  de2 <- stats::qnorm(zratio1)
  mat1 <- matrix(c(1, -r, 1/sqrt(2),
                   -r, 1, -r/sqrt(2),
                   1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  mat2 <- matrix(c(1, 0, -1/sqrt(2),
                   0, 1, -r/sqrt(2),
                   -1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res <- as.numeric(2 * (1-zratio2) * (zratio1)-
      2 * mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat1)-
      2 * mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat2))
  return(res)
}
bridgeF_tc <- function(r, zratio1, zratio2 = NULL){
  # truncated and continuous
  de1 <- stats::qnorm(zratio1)
  mat2 <- matrix(c(1, 1/sqrt(2), r/sqrt(2),
                   1/sqrt(2), 1, r,
                   r/sqrt(2), r, 1), nrow = 3)
  res <- as.numeric(-2 * fMultivar::pnorm2d(-de1, 0, rho = 1/sqrt(2)) +
                    4 * mnormt::pmnorm(c(-de1, 0, 0), mean = rep(0, 3), varcov = mat2))
  return(res)
}
bridgeF_ct <- function(r, zratio1 = NULL, zratio2){
  # continuous and truncated
  de1 <- stats::qnorm(zratio2)
  mat2 <- matrix(c(1, 1/sqrt(2), r/sqrt(2),
                   1/sqrt(2), 1, r,
                   r/sqrt(2), r, 1), nrow = 3)
  res <- as.numeric(-2 * fMultivar::pnorm2d(-de1, 0, rho = 1/sqrt(2)) +
                    4 * mnormt::pmnorm(c(-de1, 0, 0), mean = rep(0, 3), varcov = mat2))
  return(res)
}
bridgeF_bb <- function(r, zratio1, zratio2){
  # binary and binary
  de1 <- stats::qnorm(zratio1)
  de2 <- stats::qnorm(zratio2)
  res <- as.numeric(2 * (fMultivar::pnorm2d(de1, de2, rho = r) - zratio1*zratio2))
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
bridgeF_nc <- function(r, zratio1, zratio2 = NULL){
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
bridgeF_cn <- function(r, zratio1 = NULL, zratio2){
  # continuous and ternary
  de1 <- stats::qnorm(zratio2[1])
  de2 <- stats::qnorm(zratio2[2])
  mat <- matrix(c(1, 0, r/sqrt(2),
                  0, 1, -r/sqrt(2),
                  r/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res <- as.numeric(4 * fMultivar::pnorm2d(de2, 0, rho = r/sqrt(2)) - 2 * zratio2[2] +
         4 * mnormt::pmnorm(c(de1, de2, 0), mean = rep(0, 3), varcov = mat) -
         2 * zratio2[1] * zratio2[2])
  return(res)
}
bridgeF_nn <- function(r, zratio1, zratio2){
  # ternary and ternary
  de1 <- stats::qnorm(zratio1)
  de2 <- stats::qnorm(zratio2)

  res <- as.numeric(2 * fMultivar::pnorm2d(de1[2], de2[2], rho = r) * fMultivar::pnorm2d(-de1[1], -de2[1], rho = r) -
         2 * (zratio1[2] - fMultivar::pnorm2d(de1[2], de2[1], rho = r)) * (zratio2[2] - fMultivar::pnorm2d(de1[1], de2[2], rho = r)))
  return(res)
}
bridge_bn <- function(r, zratio1, zratio2){
  # binary and ternary
  de1 <- stats::qnorm(zratio1)
  de2 <- stats::qnorm(zratio2)

  res <- as.numeric(2 * fMultivar::pnorm2d(de1, de2[2], rho = r) * (1 - zratio2[1]) -
         2 * zratio2[2] * (zratio1 - fMultivar::pnorm2d(de1, de2[1], rho = r)))
  return(res)
}
bridge_nb <- function(r, zratio1, zratio2){
  # ternary and binary
  de1 <- stats::qnorm(zratio2)
  de2 <- stats::qnorm(zratio1)

  res <- as.numeric(2 * fMultivar::pnorm2d(de1, de2[2], rho = r) * (1 - zratio1[1]) -
         2 * zratio1[2] * (zratio2 - fMultivar::pnorm2d(de1, de2[1], rho = r)))
  return(res)
}
bridgeF_oc <- function(r, zratio1, zratio2 = NULL){
  # p-level ordinal and continuous
  p <- length(zratio1) + 1
  de <- stats::qnorm(zratio1)
  mat <- matrix(c(1, 0, r/sqrt(2),
                  0, 1, -r/sqrt(2),
                  r/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res = rep(NA, (p-1))
  for (i in 1:(p-1)){
    res[i] = as.numeric(4 * mnormt::pmnorm(c(de[i], de[i + 1], 0), mean = rep(0, 3), varcov = mat) -
             2 * zratio1[i] * zratio1[i + 1])
  }
  res_sum = sum(res)
  return(res_sum)
}
bridgeF_co <- function(r, zratio1 = NULL, zratio2){
  # continuous and p-level ordinal
  p <- length(zratio2) + 1
  de <- stats::qnorm(zratio2)
  mat <- matrix(c(1, 0, r/sqrt(2),
                  0, 1, -r/sqrt(2),
                  r/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res = rep(NA, (p-1))
  for (i in 1:(p-1)){
    res[i] = as.numeric(4 * mnormt::pmnorm(c(de[i], de[i + 1], 0), mean = rep(0, 3), varcov = mat) -
             2 * zratio2[i] * zratio2[i + 1])
  }
  res_sum = sum(res)
  return(res_sum)
}
