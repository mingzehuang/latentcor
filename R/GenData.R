
#' Construct a correlation matrix
#' Functions to create autocorrelation matrix (p by p) with parameter rho and block correlation matrix (p by p) using group index (of length p) and (possibly) different parameter rho for each group.
#' @rdname Sigma
#' @param p Specified matrix dimension.
#' @param rho Correlation value(s), must be between -0.99 and 0.99. Should be a scalar for \code{autocor}, and either a scalar or a vector of the same length as the maximal \code{blockind} K for \code{blockcor}.
#' @return  Correlation matrix \code{Sigma}
#' @export
autocor <- function(p, rho){
  if (abs(rho) > 0.99){ stop("correlation rho must be between -0.99 and 0.99.") }
  Sigma <- rho^abs(outer(1:p, 1:p, "-"))
  return(Sigma)
}


#' Construct a correlation matrix
#' @rdname Sigma
#' @param blockind Block index 1,\dots, K for a positive integer K specifying which variable belongs to which block, the matrix dimension is equal to \code{length(blockind)}.
#' @param rho Correlation value(s), must be between -0.99 and 0.99. Should be a scalar for \code{autocor}, and either a scalar or a vector of the same length as the maximal \code{blockind} K for \code{blockcor}.
#' @return Correlation matrix \code{Sigma}
#' @examples
#' # For p = 8,
#' # auto correlation matrix
#' autocor(8, 0.8)
#' # block correlation matrix: two blocks with the same correlation within each block
#' blockcor(c(rep(1,3), rep(2,5)), 0.8)
#' # block correlation matrix: two blocks with different correlation within each block
#' blockcor(c(rep(1,3), rep(2,5)), c(0.8, 0.3))
#'
#' @export
blockcor <- function(blockind, rho){
  if(max(abs(rho)) > 0.99){ stop("correlation rho must be between -0.99 and 0.99.") }

  p <- length(blockind)
  blk <- unique(blockind)
  if (length(rho) != length(blk)){
    if(length(rho) == 1){
      rho <- rep(rho, length(blk))
    } else {
      stop("rho and number of groups must match.")
    }
  }
  Sigma <- matrix(0, p, p)

  for (j in 1:length(blk)){
    coef <- which(blockind %in% blk[j])
    Sigma[coef, coef] <- rho[j]
  }
  diag(Sigma) = 1
  return(Sigma)
}

#' Mixed type simulation data generator
#' \code{GenData} is used to generate two sets of data of mixed types for sparse CCA under the Gaussian copula model.
#' @param n Sample size
#' @param sigma True correlation between 2 variables for bi-variable data generation.
#' @param Sigma True correlation matrix of latent variable \code{Z1} (p1 by p1).
#' @param copula1 Copula type for the first dataset. U1 = f(Z1), which could be either "exp", "cube".
#' @param copula2 Copula type for the second dataset. U2 = f(Z2), which could be either "exp", "cube".
#' @param type1 Type of the first dataset \code{X1}. Could be "continuous", "trunc", "binary", "ternary".
#' @param type2 Type of the second dataset \code{X2}. Could be "continuous", "trunc", "binary", "ternary".
#' @param muZ Mean of latent multivariate normal.
#' @param p1 Dimension of variables belong to \code{type1}.
#' @param p2 Dimension of variables belong to \code{type2}.
#' @param c1 Constant threshold for \code{X1} needed for "trunc", "binary", "ternary" and "ordinal" data type - the default is NULL.
#' @param c2 Constant threshold for \code{X2} needed for "trunc", "binary", "ternary" and "ordinal" data type - the default is NULL.
#' @return \code{GenData} returns a list containing
#' \itemize{
#'       \item{Z1: }{latent numeric data matrix (n by p1).}
#'       \item{Z2: }{latent numeric data matrix (n by p2).}
#'       \item{X1: }{observed numeric data matrix (n by p1).}
#'       \item{X2: }{observed numeric data matrix (n by p2).}
#' }
#' @export
#' @importFrom MASS mvrnorm
#' @examples
#' # Data setting
#' n <- 100; p1 <- 15; p2 <- 10 # sample size and dimensions for two datasets.
#' perm1 <- sample(1:(p1 + p2), size = p1 + p2)
#' Sigma <- autocor(p1 + p2, 0.7)[perm1, perm1]
#' mu <- rbinom(p1 + p2, 1, 0.5)
#' # Data generation
#' simdata = GenData(n = n, copula1 = "exp", copula2 = "cube", type1 = "ternary", type2 = "trunc",
#'  muZ = mu, Sigma = Sigma, p1 = p1, p2 = p2,
#'  c1 = matrix(rep(1:2, p1), nrow = 2, ncol = p1), c2 = rep(0, p2))
#'
#'
GenData <- function(n, type1, type2, sigma = NULL, p1 = 1, p2 = 1, copula1 = NULL, copula2 = NULL, muZ = NULL, Sigma = NULL, c1 = NULL, c2 = NULL){
  if (p1 == 1 & p2 == 1 & is.null(Sigma)) {
    Sigma = matrix(c(1, sigma, sigma, 1), 2, 2)
  }
  if((type1 != "continuous") & is.null(c1)){
    stop("c1 has to be defined for truncated continuous, binary, ternary or ordinal data type.")
  }
  if((type2 != "continuous") & is.null(c2)){
    stop("c2 has to be defined for truncated continuous, binary, ternary or ordinal data type.")
  }
  if (is.null(dim(c1)) & !(is.null(c1))) {
    c1 <- matrix(c1, nrow = 1, ncol = length(c1))
  }
  if (is.null(dim(c2)) & !(is.null(c2))) {
    c2 <- matrix(c2, nrow = 1, ncol = length(c2))
  }
  # jointly generate X and Y using two canonical pairs
  if (is.null(muZ)) {
    muZ <- rep(0, p1 + p2)
  }
  dat <- MASS::mvrnorm(n, mu = muZ, Sigma = Sigma) # generate a data matrix of size: n by length(muZ). length(muZ) should match with ncol(JSigma)=nrow(JSigma).

  Z1 <- as.matrix(dat[, 1:p1])
  Z2 <- as.matrix(dat[, (p1+1):(p1 + p2)])

  # Three different types of copula
    if(copula1 == "exp"){
      Z1 <- exp(Z1)
    }else if(copula1 == "cube"){
      Z1 <- Z1^3
    }

    if(copula2 == "exp"){
      Z2 <- exp(Z2)
    }else if(copula2 == "cube"){
      Z2 <- Z2^3
    }

  if(type1 != "continuous"){
    if(ncol(c1) != p1) { stop("The length of threshold vector c1 does not match with the size of the data X1.") }
    if(ncol(c1) == 1) { warning("Same threshold is applied to the all variables in the first set.") }
  }
  if(type2 != "continuous"){
    if(ncol(c2) != p2) { stop("The length of threshold vector c2 does not match with the size of the data X2.") }
    if(ncol(c2) == 1) { warning("Same threshold is applied to the all variables in the second set.") }
  }

  if(type1 == "continuous") {
    X1 <- Z1
  } else if(type1 == "trunc") {
    X1 <- ifelse(Z1 > matrix(c1, nrow = nrow(Z1), ncol = ncol(Z1), byrow = T), Z1, 0)
  } else if (type1 == "binary") {
    X1 <- ifelse(Z1 > matrix(c1, nrow = nrow(Z1), ncol = ncol(Z1), byrow = T), 1, 0)
  } else if (type1 == "ternary") {
    X1 <- ifelse(Z1 >= matrix(apply(c1, 2, max), nrow = nrow(Z1), ncol = ncol(Z1), byrow = T), 2, 1)
    X1[Z1 <= matrix(apply(c1, 2, min), nrow = nrow(Z1), ncol = ncol(Z1), byrow = T)] = 0
  }

  if(type2 == "continuous") {
    X2 <- Z2
  } else if(type2 == "trunc") {
    X2 <- ifelse(Z2 > matrix(c2, nrow = nrow(Z2), ncol = ncol(Z2), byrow = T), Z2, 0)
  } else if (type2 == "binary") {
    X2 <- ifelse(Z2 > matrix(c2, nrow = nrow(Z2), ncol = ncol(Z2), byrow = T), 1, 0)
  } else if (type2 == "ternary") {
    X2 <- ifelse(Z2 >= matrix(apply(c2, 2, max), nrow = nrow(Z2), ncol = ncol(Z2), byrow = T), 2, 1)
    X2[Z2 <= matrix(apply(c2, 2, min), nrow = nrow(Z2), ncol = ncol(Z2), byrow = T)] = 0
  }

  return(list(Z1 = Z1, Z2 = Z2, X1 = X1, X2 = X2))
}

#' Plot true correlation vs estimated correlation
#' \code{PlotPair} is to check unbiasness of estimation by plotting true correlation from simulation data vs estimated correlation from simulation data.
#' @param datapair matrix for data pairs.
#' @param namepair vector for names of data pairs.
#' @param title title for graphs.
#' @import ggplot2
#' @return \code{PlotPair} returns a plot for data1 against data2 and 45 degree benchmark line.
#' @example man/examples/estimateR_ex.R
#' @export
PlotPair <- function(datapair, namepair = c("X", "Y"), title = "Plot X vs Y") {
  df <- data.frame(datapair)
  colnames(df) = namepair
  print(ggplot(df, aes(x = datapair[ , 1], y = datapair[ , 2]))
        + geom_point(color = "blue") + geom_abline(intercept = 0, slope = 1, color = "red")
        +ggtitle(title) + xlab(namepair[1]) + ylab(namepair[2]))
}
