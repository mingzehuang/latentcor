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

GenData = function(n = 100, rhos = .5, types = c("ter", "con"), copulas = c("no", "no"), XP = NULL) {
  types = match.arg(types, c("con", "bin", "tru", "ter"), several.ok = TRUE); p = length(types)
  copulas = match.arg(copulas, c("no", "expo", "cube"), several.ok = TRUE)
  if (is.null(XP)) {
    XP = vector(mode = "list", length = p)
    for (type in unique(types)) {
      XP[types == type] = switch(type, con = NA, bin = .5, tru = .5, ter = list(c(.3, .5)))
    }
  }
  if (p == 1) {
    Z = rnorm(n = n)
    X = as.matrix(fromZtoX(z = Z, type = types, copula = copulas, xp = unlist(XP)))
  } else {
    Sigma.lower = diag(0, p); Sigma.lower[lower.tri(Sigma.lower)] = rhos
    Sigma = Sigma.lower + t(Sigma.lower) + diag(1, p)
    Z = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    X = sapply(seq(p), function(i) {fromZtoX(z = Z[ , i], type = types[i], copula = copulas[i], xp = XP[[i]])})
  }
  return(X = X)
}
