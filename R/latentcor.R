#' @title Estimate latent correlation for mixed types.
#' @description Estimation of latent correlation matrix from observed data of (possibly) mixed types (continuous/binary/truncated/ternary) based on the latent Gaussian copula model. Missing values (NA) are allowed. The estimation is based on pairwise complete observations.
#' @rdname latentcor
#' @aliases latentcor
#' @param X A numeric matrix or numeric data frame (n by p), where n is number of samples, and p is number of variables. Missing values (NA) are allowed, in which case the estimation is based on pairwise complete observations.
#' @param types A vector of length p indicating the type of each of the p variables in \code{X}. Each element must be one of \code{"con"} (continuous), \code{"bin"} (binary), \code{"ter"} (ternary) or \code{"tru"} (truncated). If the vector has length 1, then all p variables are assumed to be of the same type that is supplied. The default value is \code{NULL}, and the variable types are determined automatically using function \code{\link{get_types}}. As automatic determination of variable types takes extra time, it is recommended to supply the types explicitly when they are known in advance.
#' @param method The calculation method for latent correlations. Either \code{"original"} or \code{"approx"}. If \code{method = "approx"}, multilinear approximation method is used, which is much faster than the original method, see Yoon et al. (2021) for timing comparisons for various variable types. If \code{method = "original"}, optimization of the bridge inverse function is used. The default is \code{"approx"}.
#' @param use.nearPD Logical indicator. \code{use.nearPD = TRUE} gets nearest positive definite matrix for the estimated latent correlation matrix with shrinkage adjustment by \code{nu}. Output \code{R} is the same as \code{Rpointwise} if \code{use.nearPD = FALSE}. Default value is \code{TRUE}.
#' @param nu Shrinkage parameter for the correlation matrix, must be between 0 and 1. Guarantees that the minimal eigenvalue of returned correlation matrix is greater or equal to \code{nu}. When \code{nu = 0}, no shrinkage is performed, the returned correlation matrix will be semi-positive definite but not necessarily strictly positive definite. When \code{nu = 1}, the identity matrix is returned (not recommended).  The default (recommended) value is 0.001.
#' @param tol When \code{method = "original"}, specifies the desired accuracy of the bridge function inversion via uniroot optimization and is passed to \code{\link{optimize}}. The default value is 1e-8. When \code{method = "approx"}, this parameter is ignored.
#' @param ratio When \code{method = "approx"}, specifies the boundary value for multilinear interpolation, must be between 0 and 1. When \code{ratio = 0}, no linear interpolation is performed (the slowest execution) which is equivalent to \code{method = "original"}. When \code{ratio = 1}, linear interpolation is always performed (the fastest execution) but may lead to high approximation errors. The default (recommended) value is 0.9 which controls the approximation error and has fast execution, see Yoon et al. (2021) for details. When \code{method = "original"}, this parameter is ignored.
#' @param showplot Logical indicator. \code{showplot = TRUE} generates a ggplot object \code{plotR} with the heatmap of latent correlation matrix \code{R}. \code{plotR = NULL} if \code{showplot = FALSE}. Default value is \code{FALSE}.
#' @details
#' The function estimates latent correlation by calculating inverse bridge function (Fan et al., 2017) evaluated at the value of sample Kendall's tau (\eqn{\hat \tau}). The bridge function F connects Kendall's tau to latent correlation r so that \eqn{F(r) = E(\hat \tau)}. The form of function F depends on the variable types (continuous/binary/truncated/ternary), but is exact. The exact form of inverse is not available, so has to be evaluated numerically for each pair of variables leading to \code{Rpointwise}.
#'
#' When \code{method = "original"}, the inversion is done numerically by solving \deqn{minimize_r (F(r) - \hat \tau)^2} using \code{\link{optimize}}. The parameter \code{tol} is used to control the accuracy of the solution.
#'
#' When \code{method = "approx"}, the inversion is done approximately by interpolating previously calculated and stored values of \eqn{F^{-1}(\hat \tau)}. This is significantly faster than the original method (Yoon et al., 2021) for binary/ternary/truncated cases, however the approximation errors may be non-negligible on some regions of the space. The parameter \code{ratio} controls the region where the interpolation is performed with default recommended value of 0.9 giving a good balance of accuracy and computational speed . Increasing the value of ratio may improve speed (but possibly sacrifice the accuracy), whereas decreasing the value of ratio my improve accuracy (but possibly sacrifice the speed). See Yoon et al. 2021 and vignette for more details.
#'
#'  In case the pointwise estimator \code{Rpointwise} is has negative eigenvalues, it is projected onto the space of positive semi-definite matrices using \code{\link[Matrix]{nearPD}}. The parameter \code{nu} further allows to perform additional shrinkage towards identity matrix (desirable in cases where the number of variables p is very large) as
#'  \deqn{R = (1 - \nu) \tilde R + \nu I,}
#'  where \eqn{\tilde R} is \code{Rpointwise} after projection by \code{\link[Matrix]{nearPD}}.
#'
#' @return \code{latentcor} returns
#' \itemize{
#'       \item zratios: A list of of length p corresponding to each variable. Returns NA for continuous variable; proportion of zeros for binary/truncated variables; the cumulative proportions of zeros and ones (e.g. first value is proportion of zeros, second value is proportion of zeros and ones) for ternary variable.
#'       \item K: (p x p) Kendall Tau (Tau-a) Matrix for \code{X}
#'       \item R: (p x p) Estimated latent correlation matrix for \code{X}
#'       \item Rpointwise: (p x p) Point-wise estimates of latent correlations for \code{X}. This matrix is not guaranteed to be semi-positive definite. This is the original estimated latent correlation matrix without adjustment for positive-definiteness.
#'       \item plotR: Heatmap plot of latent correlation matrix \code{R}, NULL if \code{showplot = FALSE}
#' }
#'
#' @references
#' Fan J., Liu H., Ning Y. and Zou H. (2017) "High dimensional semiparametric latent graphical model for mixed data" \doi{10.1111/rssb.12168}.
#'
#' Yoon G., Carroll R.J. and Gaynanova I. (2020) "Sparse semiparametric canonical correlation analysis for data of mixed types" \doi{10.1093/biomet/asaa007}.
#'
#' Yoon G., Müller C.L., Gaynanova I. (2021) "Fast computation of latent correlations" \doi{10.1080/10618600.2021.1882468}.
#'
#' @import ggplot2
#' @importFrom stats quantile qnorm na.omit optimize
#' @importFrom mnormt pmnorm
#' @importFrom fMultivar pnorm2d
#' @importFrom heatmaply heatmaply
#' @importFrom Matrix nearPD
#' @importFrom pcaPP cor.fk
#' @export
#' @example man/examples/latentcor_ex.R

latentcor = function(X, types = NULL, method = c("approx", "original"), use.nearPD = TRUE, nu = 0.001, tol = 1e-8, ratio = 0.9, showplot = FALSE){
  # Check the supplied parameters are compatible with what is expected
  if(nu < 0 | nu > 1){
    stop("nu must be between 0 and 1.")
  } else if(tol <= 0) {
    stop("tol for optimization should be a positive value.")
  } else if (ratio < 0 | ratio > 1) {
    stop("ratio must be between 0 and 1.")
  }

  X = as.matrix(X)
  if (!(is.numeric(X))) {
    stop("Input data matrix should be numeric.")
  }
  p = ncol(X); name = colnames(X)

  if (is.null(types)) {
    # if variable types are not supplied, get the types by applying get_types
    types = get_types(X)
  }else{
    # if variable types are supplied, check the correct format
    types = match.arg(types, c("con", "bin", "tru", "ter"), several.ok = TRUE)
  }

  # If only one value is supplied, treat all variables as coming from the same type
  if (length(types) == 1) {
    types = rep(types, p)
  } else if (length(types) != p) {
    stop("Length of types should be either 1 for all variables or the same as number of variables (columns of X).")
  }

  method = match.arg(method, several.ok = FALSE)
  X = as.matrix(X)
  # recoding binary, truncated and ternary values.
  encodeX = encodeX(X, types)
  X = encodeX$X; types = encodeX$types; negate = encodeX$negate
  R = matrix(0, p, p); cp = rbind(row(R)[lower.tri(R)], col(R)[lower.tri(R)]); cp.col = ncol(cp)
  if (any(is.na(X))) {
    K_a.lower = sapply(seq(cp.col), function(i) Kendalltau(X[ , cp[ , i]]))
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
  R = R * outer(negate, negate)
  # Save values from pointwise estimation
  Rpointwise = R
  if (use.nearPD) {
    # Get nearest positive definite matrix.
    R = as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
    # Do adjustmnet by nu - makes it strictly positive definite like ridge
    R = (1 - nu) * R + nu * diag(nrow(R))
  }
  colnames(K) = rownames(K) = colnames(R) = rownames(R) = colnames(Rpointwise) = rownames(Rpointwise) = make.names(c(name))
  plotR = NULL
  if (showplot) {
    plotR = heatmaply(R, dendrogram = "none", main = paste0("Estimated latent correlation (", method, ")"), margins = c(80,80,80,80),
                      grid_color = "white", grid_width = 0.00001, label_names = c("Horizontal axis:", "Vertical axis:", "Latent correlation:"), limits = c(-1, 1))
  }
  return(list(zratios = zratios, K = K, R = R, Rpointwise = Rpointwise, plotR = plotR))
}
