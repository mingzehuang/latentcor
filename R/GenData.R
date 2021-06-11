
#' Construct a correlation matrix
#' Functions to create autocorrelation matrix (p by p) with parameter rho and block correlation matrix (p by p) using group index (of length p) and (possibly) different parameter rho for each group.
#' @rdname Sigma
#' @param p Specified matrix dimension.
#' @param rho Correlation value(s), must be between -0.99 and 0.99. Should be a scalar for \code{autocor}, and either a scalar or a vector of the same length as the maximal \code{blockind} K for \code{blockcor}.
#' @return  Correlation matrix \code{Sigma}
#' @export
autocor = function(p, rho){
  if (abs(rho) > 0.99){stop("correlation rho must be between -0.99 and 0.99.")}
  return(rho^abs(outer(1:p, 1:p, "-")))
}

#' Mixed type simulation data generator
#' \code{GenData} is used to generate two sets of data of mixed types for sparse CCA under the Gaussian copula model.
#' @param n Sample size
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
GenData = function(n, type1, type2, p1 = 1, p2 = 1, copula1 = NULL, copula2 = NULL, muZ = NULL, Sigma = NULL, c1 = NULL, c2 = NULL){
  if((type1 != "continuous") & is.null(c1)){
    stop("c1 has to be defined for truncated, binary, ternary data type.")
  } else if((type2 != "continuous") & is.null(c2)){
    stop("c2 has to be defined for truncated, binary, ternary data type.")
  }
  if (is.null(dim(c1)) & !(is.null(c1))) {
    c1 = matrix(c1, nrow = 1, ncol = length(c1))
  }
  if (is.null(dim(c2)) & !(is.null(c2))) {
    c2 = matrix(c2, nrow = 1, ncol = length(c2))
  }
  # jointly generate X and Y using two canonical pairs
  if (is.null(muZ)) {
    muZ = rep(0, p1 + p2)
  }
  dat = MASS::mvrnorm(n, mu = muZ, Sigma = Sigma) # generate a data matrix of size.
  Z1 = as.matrix(dat[ , 1:p1]); Z2 = as.matrix(dat[ , (p1 + 1):(p1 + p2)])
  # Three different types of copula
  X1 = fromZtoX(Z = Z1, copula = copula1, type = type1, c = c1)
  X2 = fromZtoX(Z = Z2, copula = copula2, type = type2, c = c2)
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
PlotPair = function(datapair, namepair = c("X", "Y"), title = "Plot X vs Y") {
  df = data.frame(datapair)
  colnames(df) = namepair
  print(ggplot(df, aes(x = datapair[ , 1], y = datapair[ , 2]))
        + geom_point(color = "blue") + geom_abline(intercept = 0, slope = 1, color = "red")
        +ggtitle(title) + xlab(namepair[1]) + ylab(namepair[2]))
}
