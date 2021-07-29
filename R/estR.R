#' @title Estimate latent correlation for mixed types.
#' @description Estimation of latent correlation matrix from observed data of (possibly) mixed types (continuous/binary/truncated/ternary) based on the latent Gaussian copula model. Missing values (NA) are allowed. The estimation is based on pairwise complete observations.
#' @rdname estR
#' @aliases estR
#' @param X A numeric data matrix (n by p), where n is number of samples, and p is number of variables. Missing values (NA) are allowed, in which case the estimation is based on pairwise complete observations.
#' @param types A vector of length p indicating the type of each of the p variables in \code{X}. Each element must be one of \code{"con"} (continuous), \code{"bin"} (binary), \code{"tru"} (truncated) or \code{"ter"} (ternary). If the vector has length 1, then all p variables are assumed to be of the same type that is supplied. The default value is \code{"con"} which means all variables are continuous.
#' @param method The calculation method for latent correlations. Either \code{"original"} or \code{"approx"}. If \code{method = "approx"}, multilinear approximation method is used, which is much faster than the original method. If \code{method = "original"}, optimization of the bridge inverse function is used. The default is \code{"approx"}.
#' @param nu Shrinkage parameter for the correlation matrix, must be between 0 and 1, the default value is 0.01.
#' @param tol Desired accuracy when calculating the solution of bridge function. The default value is 1e-8.
#' @param ratio The boundary value for multilinear interpolation, must be between 0 and 1, the default value is 0.9. If \code{method = "original"}, this parameter is ignored.
#' @param showplot Logical indicator. \code{showplot = TRUE} generates a ggplot object \code{plotR} with the heatmap of latent correlation matrix \code{R}. \code{plotR = NULL} if \code{showplot = FALSE}.
#' @return \code{estR} returns
#' \itemize{
#'       \item{zratios: }{A list of of length p corresponding to each variable. Returns NA for continuous variable; proportion of zeros for binary/truncated variables; the cumulative proportions of zeros and ones (e.g. first value is proportion of zeros, second value is proportion of zeros and ones) for ternary variable. }
#'       \item{K: }{(p x p) Kendall Tau (Tau-a) Matrix for \code{X} }
#'       \item{R: }{(p x p) Estimated latent correlation matrix for \code{X} }
#'       \item{Rpointwise: }{(p x p) Point-wise estimates of latent correlations for \code{X}. This matrix is not guaranteed to be semi-positive definite. If it is however, and if \code{nu=0}, then it will be the same as \code{R}.}
#'       \item{plotR: }{Heatmap plot of latent correlation matrix \code{R}, NULL if \code{showplot = FALSE}}
#' }
#'
#' @references
#' Fan J., Liu H., Ning Y. and Zou H. (2017) "High dimensional semiparametric latent graphicalmodel for mixed data" \doi{10.1111/rssb.12168}.
#'
#' Yoon G., Carroll R.J. and Gaynanova I. (2020) "Sparse semiparametric canonical correlation analysis for data of mixed types" \doi{10.1093/biomet/asaa007}.
#'
#' Yoon G., Müller C.L., Gaynanova I. (2020) "Fast computation of latent correlations" \doi{10.1080/10618600.2021.1882468}.
#'
#' @import ggplot2
#' @importFrom stats quantile qnorm na.omit optimize
#' @importFrom mnormt pmnorm
#' @importFrom fMultivar pnorm2d
#' @importFrom heatmaply heatmaply
#' @importFrom Matrix nearPD
#' @importFrom pcaPP cor.fk
#' @importFrom chebpol ipol
#' @export
#' @example man/examples/estR_ex.R

estR = function(X, types = "con", method = c("approx", "original"), nu = 0.01, tol = 1e-8, ratio = 0.9, showplot = FALSE){
  if(nu < 0 | nu > 1){
    stop("nu must be be between 0 and 1.")
  } else if(tol <= 0) {
    stop("tol for optimization should be positive value.")
  } else if (ratio < 0 | ratio > 1) {
    stop("ratio for approximation should be between 0 and 1.")
  }
  X = as.matrix(X); p = ncol(X)
  types = match.arg(types, c("con", "bin", "tru", "ter"), several.ok = TRUE)
# types = match.arg(types, c("con", "bin", "tru", "ter", "qua", "qui", "sen", "sep", "oct", "nov", "den", "dtr"), several.ok = TRUE)
  if (length(types) == 1) {
    types = rep(types, p)
  } else if (length(types) != p) {
    stop("Length of types should be either 1 or p which specifies types of p variables.")
  }
  method = match.arg(method, several.ok = FALSE)
  if (length(types) != p) {
    stop("types should have the same length as the number of variables (columns of X).")
    }
  if (length(colnames(X)) == p) {
    name = colnames(X)
  } else {
    name = paste0("X", 1:p)
  }
  R = matrix(0, p, p); cp = rbind(row(R)[lower.tri(R)], col(R)[lower.tri(R)]); cp.col = ncol(cp)
  if (any(is.na(X))) {
    K_a.lower = sapply(seq(p), function(i) Kendalltau(X[ , cp[ , i]]))
  } else {
    K_a.lower = Kendalltau(X)
  }
  zratios = zratios(X = X, types = types)
# types_code = match(types, c("con", "bin", "tru", "ter", "qua", "qui", "sen", "sep", "oct", "nov", "den", "dtr")) - 1
  types_code = match(types, c("con", "bin", "tru", "ter")) - 1
  types_cp = matrix(types_code[cp], nrow = 2); zratios_cp = matrix(zratios[cp], nrow = 2)
  types_mirror = types_cp[1, ] < types_cp[2, ]
  types_cp[ , types_mirror] = rbind(types_cp[2, types_mirror], types_cp[1, types_mirror])
  zratios_cp[ , types_mirror] = rbind(zratios_cp[2, types_mirror], zratios_cp[1, types_mirror])
  combs_cp = paste0(types_cp[1, ], types_cp[2, ]); combs = unique(combs_cp); R.lower = rep(NA, cp.col)
  for (comb in combs) {
    comb_select = combs_cp == comb
    if (comb == "00") {
      R.lower[comb_select] = sin((pi / 2) * K_a.lower[comb_select])
    } else {
      comb_select.len = sum(comb_select); K = K_a.lower[comb_select]
      zratio1 = matrix(unlist(zratios_cp[1, comb_select]), ncol = comb_select.len)
      zratio2 = matrix(unlist(zratios_cp[2, comb_select]), ncol = comb_select.len)
      R.lower[comb_select] = r_switch(method = method, K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = tol, ratio = ratio)
    }
  }
  K = matrix(0, p, p)
  K[lower.tri(K)] = K_a.lower; K = K + t(K) + diag(p); R[lower.tri(R)] = R.lower; R = R + t(R) + diag(p)
  # Save values from pointwise estimation
  Rpointwise = R
  # Check if the matrix is semi-pos.definite
  R_min_eig = min(eigen(R)$values)
  if (R_min_eig < 0) {
    message("Using Matrix::nearPD since Minimum eigenvalue of latent correlation matrix is ", R_min_eig, "smaller than 0.")
    R = as.matrix(Matrix::nearPD(R, corr = TRUE, maxit = 1000)$mat)
  }
  # Do adjustmnet by nu - makes it stricly positive definite like ridge
  R = (1 - nu) * R + nu * diag(nrow(R))
  colnames(K) = rownames(K) = colnames(R) = rownames(R) = make.names(c(name))
  plotR = NULL
  if (showplot) {
    plotR = heatmaply(R, dendrogram = "none", main = "Latent Correlation", margins = c(80,80,80,80),
                      grid_color = "white", grid_width = 0.00001, label_names = c("Horizontal axis:", "Vertical axis:", "Latent correlation:"))
  }
  return(list(zratios = zratios, K = K, R = R, Rpointwise = Rpointwise, plotR = plotR))
}
