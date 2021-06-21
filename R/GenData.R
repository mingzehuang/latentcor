#'
#' @title Mixed type simulation data generator
#' @description \code{GenData} is used to generate two sets of data of mixed types for sparse CCA under the Gaussian copula model.
#' @param n Sample size
#' @param types Types of the  dataset \code{X}. Could be "continuous", "trunc", "binary", "ternary".
#' @param rhos True correlation for auto correlated data.
#' @param copulas Copula types for the first dataset. U1 = f(Z1), which could be either "exp", "cube".
#' @param pi1 Proportions for \code{X[ , 1]} needed for "trunc", "binary", "ternary" and "ordinal" data type.
#' @param pi2 Proportions for \code{X[ , 2]} needed for "trunc", "binary", "ternary" and "ordinal" data type.
#' @return \code{GenData} returns a list containing
#' \itemize{
#'       \item{X: }{observed numeric data matrix (n by p).}
#' }
#' @export
#' @importFrom MASS mvrnorm
#' @examples
#' # Data setting
#' n = 100; p1 = 15; p2 = 10 # sample size and dimensions for two datasets.
#' # Data generation
#' GenData(n = n, type1 = "ternary", type2 = "trunc", p1 = p1, p2 = p2, rho = .9,
#'         copula1 = "cube", copula2 = "cube", c1 = c(0, 1), c2 = 0)
#'

GenData = function(n = 100, types = c("tru", "ter"), rhos = .5, copulas = NULL, pi1 = .5, pi2 = c(.3, .5), ...) {
  p = length(types); Sigma.lower = diag(0, p); Sigma.lower[lower.tri(Sigma.lower)] = rhos
  Z = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma.lower + t(Sigma.lower) + diag(1, p))
  X = sapply(seq(p), function(i) {fromZtoX(z = Z[ , i], type = types[i], copula = copulas[i], get(paste0("pi", i)))})
  return(X = X)
}

fromZtoX = function(z, type, copula, pi) {
  if (is.null(copula)) {
    u = z
  } else if(copula == "exp"){
    u = exp(z)
  }else if(copula == "cube"){
    u = z^3
  }
  if(type == "continuous") {
    x = u
  } else {
    if (is.null(pi)) {
      stop("Proportions need to be set for binary, truncated and ternary cases.")
    } else {
      q = quantile(u, cumsum(pi))
      if (type == "binary") {
        x = ifelse(u > q, 1, 0)
      } else if(type == "trunc") {
        x = ifelse(u > q, u, q) - q
      } else if (type == "ternary") {
        x = rep(1, length(u))
        x[u > q[2]] = 2; x[u <= q[1]] = 0
      }
    }
  }
  return(x)
}
