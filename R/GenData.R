#'
#' @title Mixed type simulation data generator
#' @description \code{GenData} is used to generate two sets of data of mixed types for sparse CCA under the Gaussian copula model.
#' @param n Sample size
#' @param type1 Type of the first dataset \code{X1}. Could be "continuous", "trunc", "binary", "ternary".
#' @param type2 Type of the second dataset \code{X2}. Could be "continuous", "trunc", "binary", "ternary".
#' @param p1 Dimension of variables belong to \code{type1}.
#' @param p2 Dimension of variables belong to \code{type2}.
#' @param mu Mean of latent multivariate normal.
#' @param Sigma True correlation matrix of latent variable \code{Z1} (p1 by p1).
#' @param Z1 continuous pre-generated data matrix to be converted into other data type.
#' @param Z2 continuous pre-generated data matrix to be converted into other data type.
#' @param rho True correlation for auto correlated data.
#' @param bdim Dimensions for subdiagonal matrices.
#' @param copula1 Copula type for the first dataset. U1 = f(Z1), which could be either "exp", "cube".
#' @param copula2 Copula type for the second dataset. U2 = f(Z2), which could be either "exp", "cube".
#' @param c1 Constant threshold for \code{X1} needed for "trunc", "binary", "ternary" and "ordinal" data type - the default is NULL.
#' @param c2 Constant threshold for \code{X2} needed for "trunc", "binary", "ternary" and "ordinal" data type - the default is NULL.
#' @param q1 quantile threshold for \code{X1} needed for "trunc", "binary", "ternary" and "ordinal" data type - the default is NULL.
#' @param q2 quantile threshold for \code{X2} needed for "trunc", "binary", "ternary" and "ordinal" data type - the default is NULL.
#' @return \code{GenData} returns a list containing
#' \itemize{
#'       \item{Z1: }{latent numeric data matrix (n by p1).}
#'       \item{Z2: }{latent numeric data matrix (n by p2).}
#'       \item{X1: }{observed numeric data matrix (n by p1).}
#'       \item{X2: }{observed numeric data matrix (n by p2).}
#' }
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom Matrix .bdiag
#' @examples
#' # Data setting
#' n = 100; p1 = 15; p2 = 10 # sample size and dimensions for two datasets.
#' # Data generation
#' GenData(n = n, type1 = "ternary", type2 = "trunc", p1 = p1, p2 = p2, rho = .9,
#'         copula1 = "cube", copula2 = "cube", c1 = c(0, 1), c2 = 0)
#'
#'
GenData = function(n, type1 = "continuous", type2 = NULL, p1 = 1, p2 = NULL, rho = NULL, bdim = NULL, copula1 = NULL, copula2 = NULL, mu = NULL, Sigma = NULL, Z1 = NULL, Z2 = NULL, c1 = NULL, c2 = NULL, q1 = NULL, q2 = NULL){
  if (!(is.null(Z1)) & !(is.null(Z2))) {
    X1 = fromZtoX(Z = Z1, copula = copula1, type = type1, q = q1, c = c1)
    X2 = fromZtoX(Z = Z2, copula = copula2, type = type2, q = q2, c = c2)
    return(list(X1 = X1, X2 = X2))
  } else if (!(is.null(Z1))) {
    X1 = fromZtoX(Z = Z1, copula = copula1, type = type1, q = q1, c = c1)
    return(list(X1 = X1))
  } else if (is.null(type2)) {
    if (is.null(mu)) {mu = rep(0, p1)}
    if (p1 == 1) {
      Sigma = ifelse(is.null(Sigma), 1, Sigma)
      Z1 = rnorm(n = n, mean = mu, sd = sqrt(Sigma))
    } else {
      if (is.null(Sigma) & is.null(rho)) {
        stop("correlation rho or covariance matrix Sigma needs to be specified.")
      } else if (is.null(Sigma) & is.null(bdim)) {
        Sigma = autocor(p1, rho)
      } else if (is.null(Sigma)) {
        Sigma = as.matrix(.bdiag(mapply(function(x, y){autocor(x, y)}, bdim, rho)))
      }
      Z1 = MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
    }
    X1 = fromZtoX(Z = Z1, copula = copula1, type = type1, q = q1, c = c1)
    return(list(mu = mu, Sigma = Sigma, Z1 = Z1, X1 = X1))
  } else {
    if (is.null(mu)) {mu = rep(0, p1 + p2)}
    if (is.null(Sigma) & is.null(rho)) {
      stop("correlation rho or covariance matrix Sigma needs to be specified.")
    } else if (is.null(Sigma) & is.null(bdim)) {
      Sigma = autocor(p1 + p2, rho)
    } else if (is.null(Sigma)) {
      Sigma = as.matrix(.bdiag(mapply(function(x, y){autocor(x, y)}, bdim, rho)))
    }
    dat = MASS::mvrnorm(n, mu = mu, Sigma = Sigma) # generate a data matrix of size.
    Z1 = as.matrix(dat[ , 1:p1]); Z2 = as.matrix(dat[ , (p1 + 1):(p1 + p2)])
    X1 = fromZtoX(Z = Z1, copula = copula1, type = type1, q = q1, c = c1)
    X2 = fromZtoX(Z = Z2, copula = copula2, type = type2, q = q2, c = c2)
    return(list(mu = mu, Sigma = Sigma, Z1 = Z1, Z2 = Z2, X1 = X1, X2 = X2))
  }
}

