#'
#' @title Mixed type simulation data generator
#' @description \code{GenData} is used to generate two sets of data of mixed types for sparse CCA under the Gaussian copula model.
#' @param n Sample size
#' @param types Types of the  dataset \code{X}. Could be "con" (continuous), "bin" (binary), "tru" (truncated) or "ter" (ternary).
#' @param rhos True correlation for auto correlated data.
#' @param copulas Copula types for the first dataset. U = f(Z), which could be "NA", "exp" or "cube".
#' @param XP list of proportions corresponding to each data series.
#' @return \code{GenData} returns a list containing
#' \itemize{
#'       \item{X: }{observed numeric data matrix (n by p).}
#' }
#' @export
#' @importFrom stats rnorm
#' @importFrom MASS mvrnorm
#' @examples
#' GenData()

GenData = function(n = 100, types = c("tru", "ter"), rhos = .5, copulas = c(NA, NA), XP = list(.5, c(.3, .5))) {
  p = length(types); types_code = rep(NA, p); types_code = as.numeric(type_list[types])
  if (p == 1) {
    Z = rnorm(n = n)
    X = as.matrix(fromZtoX(z = Z, type_code = types_code, copula = copulas, xp = unlist(XP)))
  } else {
    Sigma.lower = diag(0, p); Sigma.lower[lower.tri(Sigma.lower)] = rhos
    Sigma = Sigma.lower + t(Sigma.lower) + diag(1, p)
    Z = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    X = sapply(seq(p), function(i) {fromZtoX(z = Z[ , i], type_code = types_code[i], copula = copulas[i], xp = XP[[i]])})
  }
  return(X = X)
}

fromZtoX = function(z, type_code, copula, xp) {
  x = z
  if (!(is.na(copula))) {x = copula_list[[copula]](z)}
  if(type_code != 0) {x = transform_list[[type_code]](u = x, xp = xp)}
  return(x)
}
