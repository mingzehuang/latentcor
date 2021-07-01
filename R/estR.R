#' @title Calculate Kendall's Tau matrix then estimate latent correlation matrix
#' @description Estimation of latent correlation matrix from observed data of (possibly) mixed types (continuous/biary/truncated continuous) based on the latent Gaussian copula model.
#' @rdname estR
#' @aliases estR
#' @param X A numeric data matrix (n by p).
#' @param types A vector with length of p which specifies types of variables in \code{X} correspondingly.
#'              Must be one of "con" (continuous), "bin" (binary), "tru" (truncated) or "ter" (ternary). The default value is c("ter", "con") which means the first variable is ternary, second variable is continuous.
#' @param method The calculation method of latent correlation. Either "original" method or "approx". If \code{method = "approx"}, multilinear approximation method is used, which is much faster than the original method.
#'               If \code{method = "original"}, optimization of the bridge inverse function is used. The default is "approx".
#' @param nu Shrinkage parameter for correlation matrix, must be between 0 and 1, the default value is 0.01.
#' @param tol Desired accuracy when calculating the solution of bridge function. The default value is 1e-8.
#' @param ratio The maximum ratio of Kendall's tau and boundary to implement multilinear interpolation. The default value is 0.9.
#' @param showplot Plot latent correlation matrix \code{R} as a heatmap.
#' @return \code{estR} returns
#' \itemize{
#'       \item{zratios: }{A list of zratios. Each element correponds to the zratios for one variable.}
#'       \item{K: }{Kendall Tau (Tau-a) Matrix of \code{X} (p x p)}
#'       \item{R: }{Estimated latent correlation matrix of whole \code{X} (p x p)}
#'       \item{plotR: }{Heatmap plot for latent correlation matrix \code{R}}
#' }
#' @references
#' Fan J., Liu H., Ning Y. and Zou H. (2017) "High dimensional semiparametric latent graphicalmodel for mixed data" <doi:10.1111/rssb.12168>.
#' Yoon G., Carroll R.J. and Gaynanova I. (2020) "Sparse semiparametric canonical correlation analysis for data of mixed types" <doi:10.1093/biomet/asaa007>.
#' Yoon G., Müller C.L., Gaynanova I. (2020) "Fast computation of latent correlations" <arXiv:2006.13875>.
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

estR = function(X, types = c("ter", "con"), method = "approx", nu = 0.01, tol = 1e-8, ratio = 0.9, showplot = FALSE){
  if(nu < 0 | nu > 1){
    stop("nu must be be between 0 and 1.")
  } else if(tol <= 0) {
    stop("tol for optimization should be positive value.")
  } else if (ratio < 0 | ratio > 1) {
    stop("ratio for approximation should be between 0 and 1.")
  }
  types = match.arg(types, c("con", "bin", "tru", "ter", "qua", "qui", "sen", "sep", "oct", "nov", "den"), several.ok = TRUE)
  method = match.arg(method, c("original", "approx"), several.ok = FALSE)
  X = as.matrix(X); p = ncol(X);
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
  types_code = match(types, c("con", "bin", "tru", "ter", "qua", "qui", "sen", "sep", "oct", "nov", "den")) - 1
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
  R_min_eig = min(eigen(R)$values)
  if (R_min_eig < 0) {
    message("Use Matrix::nearPD since Minimum eigenvalue of latent correlation matrix is ", R_min_eig, "smaller than 0.")
    R = as.matrix(Matrix::nearPD(R, corr = TRUE, maxit = 1000)$mat)
  }
  R = (1 - nu) * R + nu * diag(nrow(R))
  colnames(K) = rownames(K) = colnames(R) = rownames(R) = make.names(c(name))
  plotR = NULL
  if (showplot) {
    plotR = heatmaply(R, dendrogram = "none", main = "Latent Correlation", margins = c(80,80,80,80),
                      grid_color = "white", grid_width = 0.00001, label_names = c("Horizontal axis:", "Vertical axis:", "Latent correlation:"))
  }
  return(list(zratios = zratios, K = K, R = R, plotR = plotR))
}
