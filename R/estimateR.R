
#' @title Estimate latent correlation matrix
#'
#' @description Estimation of latent correlation matrix from observed data of (possibly) mixed types (continuous/biary/truncated continuous) based on the latent Gaussian copula model.
#'
#' @aliases estimateR estimateR_mixed
#' @param X A numeric data matrix (n by p), n is the sample size and p is the number of variables.
#' @param type A type of variables in \code{X}, must be one of "continuous", "binary" or "trunc".
#' @param method The calculation method of latent correlation. Either "original" method or "approx". If \code{method = "approx"}, multilinear approximation method is used, which is much faster than the original method. If \code{method = "original"}, optimization of the bridge inverse function is used. The default is "approx".
#' @param use.nearPD A logical value indicating whether to use \link[Matrix]{nearPD} or not when the resulting correlation estimator is not positive definite (have at least one negative eigenvalue).
#' @param nu Shrinkage parameter for correlation matrix, must be between 0 and 1, the default value is 0.01.
#' @param tol Desired accuracy when calculating the solution of bridge function.
#' @param verbose If \code{verbose = FALSE}, printing information whether nearPD is used or not is disabled. The default value is FALSE.
#' @param ratio The maximum ratio of Kendall's tau and boundary to implement multilinear interpolation.
#' @return \code{estimateR} returns
#' \itemize{
#'       \item{type: }{Type of the data matrix \code{X}}
#'       \item{R: }{Estimated p by p latent correlation matrix of \code{X}}
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

estimateR <- function(X, type, method, use.nearPD = TRUE, nu = 0.01, tol = 1e-6, verbose = FALSE, ratio = 0.9){
  X <- as.matrix(X)
  p <- ncol(X)
  zratio = zratio(X = X, type = type)
  K = Kendall(X = X, type = type)
  R = fromKtoR(K = K, zratio = zratio, type = type, method = method, tol = tol, ratio = ratio)
  R.final = R_adj(R = R, use.nearPD = use.nearPD, verbose = verbose, nu = nu)
  ### To keep the correct column names for each matrices
  if(length(colnames(X)) == p){
    colnames(R.final) <- rownames(R.final) <- c(colnames(X))
  }
  return(list(type = type, R = R.final))
}

#'
#'
#'
#' @title Estimate latent correlation matrix
#' @rdname estimateR
#' @aliases estimateR estimateR_mixed
#' @param X1 A numeric data matrix (n by p1).
#' @param X2 A numeric data matrix (n by p2).
#' @param type1 A type of variables in \code{X1}, must be one of "continuous", "binary" or "trunc".
#' @param type2 A type of variables in \code{X2}, must be one of "continuous", "binary" or "trunc".
#' @param method The calculation method of latent correlation. Either "original" method or "approx". If \code{method = "approx"}, multilinear approximation method is used, which is much faster than the original method. If \code{method = "original"}, optimization of the bridge inverse function is used. The default is "approx".
#' @param use.nearPD A logical value indicating whether to use \link[Matrix]{nearPD} or not when the resulting correlation estimator is not positive definite (have at least one negative eigenvalue).
#' @param nu Shrinkage parameter for correlation matrix, must be between 0 and 1, the default value is 0.01.
#' @param tol Desired accuracy when calculating the solution of bridge function.
#' @param verbose If \code{verbose = FALSE}, printing information whether nearPD is used or not is disabled. The default value is FALSE.
#'
#'
#' @return \code{estimateR_mixed} returns
#' \itemize{
#'       \item{type1: }{Type of the data matrix \code{X1}}
#'       \item{type2: }{Type of the data matrix \code{X2}}
#'       \item{R: }{Estimated latent correlation matrix of whole \code{X} = (\code{X1}, \code{X2}) (p1+p2 by p1+p2)}
#'       \item{R1: }{Estimated latent correlation matrix of \code{X1} (p1 by p1)}
#'       \item{R2: }{Estimated latent correlation matrix of \code{X2} (p2 by p2)}
#'       \item{R12: }{Estimated latent correlation matrix between \code{X1} and \code{X2} (p1 by p2)}
#' }
#'
#' @export
estimateR_mixed <- function(X1, X2 = NULL, type1, type2 = NULL, method = "approx", use.nearPD = TRUE, nu = 0.01, tol = 1e-6, verbose = FALSE, ratio = 0.9){
  # shrinkage method
  if(nu < 0 | nu > 1){
    stop("nu must be be between 0 and 1.")
  }
  if (is.null(X2)) {
    out = estimateR(X = X1, type = type1, method = method, use.nearPD = use.nearPD, nu = nu, tol = tol, verbose = verbose, ratio = ratio)
    return(out)
  } else {
  X1 <- as.matrix(X1)
  p1 <- ncol(X1)
  zratio1 = zratio(X = X1, type = type1)
  X2 <- as.matrix(X2)
  p2 <- ncol(X2)
  if (nrow(X1) != nrow(X2)){ # Check of they have the same sample size.
    stop ("X1 and X2 must have the same sample size.")
  }
  zratio2 = zratio(X = X2, type = type2)

  if (p1 == 1 & p2 == 1){
    # This is just pairwise correlation
    k12 = KendallTau(X1, X2)
    r12 = fromKtoR_mixed(K12 = k12, zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2, method = method, tol = tol, ratio = ratio)
    return(list(type = c(type1, type2), R1 = 1, R2 = 1, R12 = r12, R = matrix(c(1, r12, r12, 1), 2, 2)))
  }

  if (type1 == type2) {
    ################### both datasets are of the same type. CC, TT or BB case.
    Xcomb <- cbind(X1, X2)
    R.final <- estimateR(Xcomb, type = type1, method = method, use.nearPD = use.nearPD, nu = nu, tol = tol, ratio = ratio)$R
    R1 <- R.final[1:p1, 1:p1]
    R2 <- R.final[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
    R12 <- R.final[1:p1, (p1 + 1):(p1 + p2)]
  } else {
    ################### datasets are of different type
    # Start with 1st dataset
    if (p1 == 1){
      R1 <- as.matrix(1)
    }else {
      K1 = Kendall(X = X1, type = type1)
      R1 = fromKtoR(K = K1, zratio = zratio1, type = type1, method = method, tol = tol, ratio = ratio)
    }
    # Continue with 2nd dataset
    if (p2 == 1){
      R2 <- as.matrix(1)
    }else {
      K2 = Kendall(X = X2, type = type2)
      R2 = fromKtoR(K = K2, zratio = zratio2, type = type2, method = method, tol = tol, ratio = ratio)
    }
    # Do cross-product
    K12 <- Kendall_matrix(X1, X2)
    R12 <- fromKtoR_mixed(K12, zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2, method = method, tol = tol, ratio = ratio)

    Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))

    R.final = R_adj(R = Rall, use.nearPD = use.nearPD, verbose = verbose, nu = nu)

    ### To keep the column names in R according to column names that are originally supplied in each matrix
    if(length(colnames(X1)) == p1 & length(colnames(X2)) == p2){
      colnames(R.final) <- rownames(R.final) <- c(colnames(X1), colnames(X2))
    } else if(length(colnames(X1)) != p1 & length(colnames(X2)) == p2){
      colnames(R.final) <- rownames(R.final) <- rep(NA, p1 + p2)
      colnames(R.final)[(p1 + 1):(p1 + p2)] <- rownames(R.final)[(p1 + 1):(p1 + p2)] <- colnames(X2)
    } else if(length(colnames(X1)) == p1 & length(colnames(X2)) != p2){
      colnames(R.final) <- rownames(R.final) <- rep(NA, p1 + p2)
      colnames(R.final)[1:p1] <- rownames(R.final)[1:p1] <- colnames(X1)
    }

    # For convenience, split the R matrices
    R1 <- R.final[1:p1, 1:p1]
    R2 <- R.final[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
    R12 <- R.final[1:p1, (p1 + 1):(p1 + p2)]
  }
  return(list(type = c(type1, type2), R1 = R1, R2 = R2, R12 = R12, R = R.final))
  }
}
