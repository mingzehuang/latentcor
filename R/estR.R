#'
#'
#'
#' @title Estimate latent correlation matrix
#' @description Estimation of latent correlation matrix from observed data of (possibly) mixed types (continuous/biary/truncated continuous) based on the latent Gaussian copula model.
#' @rdname estR
#' @aliases estR
#' @param X1 A numeric data matrix (n by p1).
#' @param X2 A numeric data matrix (n by p2).
#' @param type1 A type of variables in \code{X1}, must be one of "continuous", "binary" or "trunc".
#' @param type2 A type of variables in \code{X2}, must be one of "continuous", "binary" or "trunc".
#' @param method The calculation method of latent correlation. Either "original" method or "approx". If \code{method = "approx"}, multilinear approximation method is used, which is much faster than the original method. If \code{method = "original"}, optimization of the bridge inverse function is used. The default is "approx".
#' @param use.nearPD A logical value indicating whether to use \link[Matrix]{nearPD} or not when the resulting correlation estimator is not positive definite (have at least one negative eigenvalue).
#' @param nu Shrinkage parameter for correlation matrix, must be between 0 and 1, the default value is 0.01.
#' @param tol Desired accuracy when calculating the solution of bridge function.
#' @param verbose If \code{verbose = FALSE}, printing information whether nearPD is used or not is disabled. The default value is FALSE.
#' @param ratio The maximum ratio of Kendall's tau and boundary to implement multilinear interpolation.
#'
#' @return \code{estR} returns
#' \itemize{
#'       \item{type1: }{Type of the data matrix \code{X1}}
#'       \item{type2: }{Type of the data matrix \code{X2}}
#'       \item{R: }{Estimated latent correlation matrix of whole \code{X} = (\code{X1}, \code{X2}) (p1+p2 by p1+p2)}
#'       \item{R1: }{Estimated latent correlation matrix of \code{X1} (p1 by p1)}
#'       \item{R2: }{Estimated latent correlation matrix of \code{X2} (p2 by p2)}
#'       \item{R12: }{Estimated latent correlation matrix between \code{X1} and \code{X2} (p1 by p2)}
#' }
#' @references
#' Fan J., Liu H., Ning Y. and Zou H. (2017) "High dimensional semiparametric latent graphicalmodel for mixed data" <doi:10.1111/rssb.12168>.
#'
#' Yoon G., Carroll R.J. and Gaynanova I. (2020) "Sparse semiparametric canonical correlation analysis for data of mixed types" <doi:10.1093/biomet/asaa007>.
#'
#' Yoon G., Müller C.L., Gaynanova I. (2020) "Fast computation of latent correlations" <arXiv:2006.13875>.
#'
#' @export
#' @example man/examples/estimateR_ex.R
#'

estR = function(X1, type1, X2 = NULL, type2 = NULL, method = "approx", use.nearPD = TRUE, nu = 0.01, tol = 1e-8, verbose = FALSE, ratio = 0.9){
  # shrinkage method
  if(nu < 0 | nu > 1){
    stop("nu must be be between 0 and 1.")
  }
  X1 = as.matrix(X1); p1 = ncol(X1)
  if (length(colnames(X1)) == p1) {
    name1 = colnames(X1)
  } else {
    name1 = paste("X1", 1:p1, sep = "")
  }
  if (is.null(X2)) {
    if (p1 == 1) {
      R = as.matrix(1)
    } else {
      R = estimateR(X1 = X1, type1 = type1, method = method, tol = tol, ratio = ratio)
    }
    colnames(R) = rownames(R) = make.names(c(name1))
    return(R)
  } else {
    X2 = as.matrix(X2); p2 = ncol(X2)
    if (length(colnames(X2)) == p2) {
      name2 = colnames(X2)
    } else {
      name2 = paste("X2", 1:p2, sep = "")
    }
    if (nrow(X1) != nrow(X2)){ # Check of they have the same sample size.
    stop ("X1 and X2 must have the same sample size.")
    } else {
      R1 = estimateR(X1 = X1, type1 = type1, method = method, tol = tol, ratio = ratio)
      R2 = estimateR(X1 = X2, type1 = type2, method = method, tol = tol, ratio = ratio)
      R12 = estimateR(X1 = X1, type1 = type1, X2 = X2, type2 = type2, method = method, tol = tol, ratio = ratio)
      Rall = rbind(cbind(R1, R12), cbind(t(R12), R2))
      R.final = R_adj(R = Rall, use.nearPD = use.nearPD, verbose = verbose, nu = nu)
    }
    ### To keep the column names in R according to column names that are originally supplied in each matrix
    colnames(R.final) = rownames(R.final) = make.names(c(name1, name2))
    # For convenience, split the R matrices
    R1 = R.final[1:p1, 1:p1]
    R2 = R.final[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
    R12 = R.final[1:p1, (p1 + 1):(p1 + p2)]
    return(list(type = c(type1, type2), R1 = R1, R2 = R2, R12 = R12, R = R.final))
  }
}
