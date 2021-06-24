#' @title Estimate latent correlation matrix
#' @description Estimation of latent correlation matrix from observed data of (possibly) mixed types (continuous/biary/truncated continuous) based on the latent Gaussian copula model.
#' @rdname estR
#' @aliases estR
#' @param X A numeric data matrix (n by p).
#' @param types A type of variables in \code{X}, must be one of "continuous", "binary", "trunc" or "ternary".
#' @param method The calculation method of latent correlation. Either "original" method or "approx". If \code{method = "approx"}, multilinear approximation method is used, which is much faster than the original method. If \code{method = "original"}, optimization of the bridge inverse function is used. The default is "approx".
#' @param nu Shrinkage parameter for correlation matrix, must be between 0 and 1, the default value is 0.01.
#' @param tol Desired accuracy when calculating the solution of bridge function.
#' @param ratio The maximum ratio of Kendall's tau and boundary to implement multilinear interpolation.
#' @param corplot Plot latent correlation matrix \code{R} as a heatmap.
#' @return \code{estR} returns
#' \itemize{
#'       \item{K: }{Kendall Tau Matrix of \code{X} (p x p)}
#'       \item{R: }{Estimated latent correlation matrix of whole \code{X} (p x p)}
#' }
#' @references
#' Fan J., Liu H., Ning Y. and Zou H. (2017) "High dimensional semiparametric latent graphicalmodel for mixed data" <doi:10.1111/rssb.12168>.
#' Yoon G., Carroll R.J. and Gaynanova I. (2020) "Sparse semiparametric canonical correlation analysis for data of mixed types" <doi:10.1093/biomet/asaa007>.
#' Yoon G., Müller C.L., Gaynanova I. (2020) "Fast computation of latent correlations" <arXiv:2006.13875>.
#' @import ggplot2
#' @importFrom utils combn
#' @importFrom stats qnorm na.omit optimize
#' @importFrom mnormt pmnorm
#' @importFrom fMultivar pnorm2d
#' @importFrom heatmaply heatmaply
#' @importFrom Matrix nearPD
#' @importFrom pcaPP cor.fk
#' @importFrom chebpol ipol
#' @export
#' @example man/examples/estR_ex.R

estR = function(X, types = c("ter", "con"), method = "approx", nu = 0.01, tol = 1e-8, ratio = 0.9, corplot = FALSE){
  # shrinkage method
  if(nu < 0 | nu > 1){
    stop("nu must be be between 0 and 1.")
  }
  types = match.arg(types, c("con", "bin", "tru", "ter"), several.ok = TRUE)
  method = match.arg(method, c("original", "approx"), several.ok = FALSE)
  X = as.matrix(X); X = na.omit(X); n = nrow(X); p = ncol(X);
  if (length(colnames(X)) == p) {
    name = colnames(X)
  } else {
    name = paste0("X", 1:p)
  }
  types_code = as.numeric(type_list[types])
  cp = combn(p, 2); cp.col = ncol(cp); n0 = n * (n - 1) / 2
  n_X = apply(X, 2, function(x) {n_x(x = x, n)})
  n_X_sqd = sqrt(n0 - n_X)
  btoa_cp = matrix(n_X_sqd[cp], nrow = 2)
  btoa2 = btoa_cp[1, ] * btoa_cp[2, ] / n0
  K_b = pcaPP::cor.fk(X)
  K_b.lower = K_b[lower.tri(K_b)]
  K_a.lower = K_b.lower * btoa2
  types_cp = matrix(types_code[cp], nrow = 2)
  zratios = zratios(X = X, types_code = types_code)
  zratios_cp = matrix(zratios[cp], nrow = 2)
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
      zratio2 = unlist(zratios_cp[2, comb_select])
      if (!(is.null(zratio2))) {zratio2 = matrix(zratio2, ncol = comb_select.len)}
      if (method == "original") {
        R.lower[comb_select] = r_sol(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = tol)
      } else {
        cutoff = abs(K) > ratio * bound_list[[comb]](zratio1 = zratio1, zratio2 = zratio2)
        if (sum(cutoff) == 0) {
          R.lower[comb_select] = r_ml(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb)
        } else if (sum(!(cutoff)) == 0) {
          R.lower[comb_select] = r_sol(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = tol)
        } else {
          R.lower[comb_select][cutoff] = r_sol(K = K[cutoff], zratio1 = zratio1[ , cutoff], zratio2 = zratio2[ , cutoff], comb = comb, tol = tol)
          R.lower[comb_select][!(cutoff)] = r_ml(K = K[!(cutoff)], zratio1 = zratio1[ , !(cutoff)], zratio2 = zratio2[ , !(cutoff)], comb = comb)
        }
      }
    }
  }
  K = R = matrix(0, p, p)
  K[lower.tri(K)] = K_a.lower; K = K + t(K) + diag(p); R[lower.tri(R)] = R.lower; R = R + t(R) + diag(p)
  R_min_eig = min(eigen(R)$values)
  if (R_min_eig < 0) {
    message("Use Matrix::nearPD since Minimum eigenvalue of latent correlation matrix is ", R_min_eig, "smaller than 0.")
    R = as.matrix(Matrix::nearPD(R, corr = TRUE, maxit = 1000)$mat)
  }
  R = (1 - nu) * R + nu * diag(nrow(R))
  colnames(K) = rownames(K) = colnames(R) = rownames(R) = make.names(c(name))
  plotR = NULL
  if (corplot) {
    plotR = heatmaply(R, dendrogram = "none", main = "Latent Correlation", margins = c(80,80,80,80),
                      grid_color = "white", grid_width = 0.00001, label_names = c("Horizontal axis:", "Vertical axis:", "Latent correlation:"))
  }
  return(list(K = K, R = R, plotR = plotR))
}

n_x = function(x, n) {
  if (length(unique(x) != n)) {
    x.info = rle(sort(x))
    t_x = x.info$lengths[x.info$lengths > 1]
    n_x = sum(t_x * (t_x - 1) / 2)
  } else {
    n_x = 0
  }
  return(n_x)
}

zratios = function(X, types_code) {
  X = as.matrix(X)
  out = vector(mode = "list", length = ncol(X))
  for (i in unique(types_code[types_code != 0])) {
    out[types_code == i] = zratio_list[[i]](X[ , types_code == i])
  }
  return(out)
}

r_sol = function(K, zratio1, zratio2, comb, tol) {
  bridge_comb = bridge_list[[comb]]
  K.len = length(K); out = rep(NA, K.len)
  zratio1 = as.matrix(zratio1)
  if (!(is.null(zratio2))) {zratio2 = as.matrix(zratio2)}
  for (i in K.len) {
    f = function(r)(bridge_comb(r = r, zratio1 = zratio1[ , i], zratio2 = zratio2[ , i]) - K[i])^2
    op = tryCatch(optimize(f, lower = -0.999, upper = 0.999, tol = tol)[1], error = function(e) 100)
    if(op == 100) {
      warning("Optimize returned error one of the pairwise correlations, returning NA")
      out[i] = NA
    } else {
    out[i] = unlist(op)
    }
  }
  return(out)
}

r_ml = function(K, zratio1, zratio2, comb) {
  ipol_comb = ipol_list[[comb]]
  zratio1 = as.matrix(zratio1); zratio1.row = nrow(zratio1)
  if (!(is.null(zratio2))) {
    zratio2 = as.matrix(zratio2); zratio2.row = nrow(zratio2)
  } else {
    zratio2.row = 0
  }
  K = K / bound_list[[comb]](zratio1 = zratio1, zratio2 = zratio2)
  if (zratio1.row > 1) {zratio1[1:(zratio1.row - 1), ] = zratio1[1:(zratio1.row - 1), ] / zratio1[2:zratio1.row, ]}
  if (zratio2.row > 1) {zratio2[1:(zratio2.row - 1), ] = zratio2[1:(zratio2.row - 1), ] / zratio2[2:zratio2.row, ]}
  out = ipol_comb(rbind(K, zratio1, zratio2)) / 10^7
  return(out)
}
