#' @title Estimate latent correlation matrix
#' @description Estimation of latent correlation matrix from observed data of (possibly) mixed types (continuous/biary/truncated continuous) based on the latent Gaussian copula model.
#' @rdname estR
#' @aliases estR
#' @param X A numeric data matrix (n by p).
#' @param type A type of variables in \code{X}, must be one of "continuous", "binary", "trunc" or "ternary".
#' @param method The calculation method of latent correlation. Either "original" method or "approx". If \code{method = "approx"}, multilinear approximation method is used, which is much faster than the original method. If \code{method = "original"}, optimization of the bridge inverse function is used. The default is "approx".
#' @param use.nearPD A logical value indicating whether to use \link[Matrix]{nearPD} or not when the resulting correlation estimator is not positive definite (have at least one negative eigenvalue).
#' @param nu Shrinkage parameter for correlation matrix, must be between 0 and 1, the default value is 0.01.
#' @param tol Desired accuracy when calculating the solution of bridge function.
#' @param verbose If \code{verbose = FALSE}, printing information whether nearPD is used or not is disabled. The default value is FALSE.
#' @param ratio The maximum ratio of Kendall's tau and boundary to implement multilinear interpolation.
#' @param plotR Plot latent correlation matrix \code{R} as a heatmap.
#' @return \code{estR} returns
#' \itemize{
#'       \item{R: }{Estimated latent correlation matrix of whole \code{X} (p by p)}
#' }
#' @references
#' Fan J., Liu H., Ning Y. and Zou H. (2017) "High dimensional semiparametric latent graphicalmodel for mixed data" <doi:10.1111/rssb.12168>.
#' Yoon G., Carroll R.J. and Gaynanova I. (2020) "Sparse semiparametric canonical correlation analysis for data of mixed types" <doi:10.1093/biomet/asaa007>.
#' Yoon G., Müller C.L., Gaynanova I. (2020) "Fast computation of latent correlations" <arXiv:2006.13875>.
#' @import stats ggplot2
#' @importFrom heatmaply heatmaply
#' @importFrom Matrix nearPD
#' @importFrom pcaPP cor.fk
#' @importFrom chebpol ipol
#' @export
#' @example man/examples/estR_ex.R
estR = function(X, types, method = "approx", nu = 0.01, tol = 1e-8, ratio = 0.9, plotR = FALSE){
  # shrinkage method
  if(nu < 0 | nu > 1){
    stop("nu must be be between 0 and 1.")
  }
  X = as.matrix(X); X = na.omit(X); n = nrow(X); p = ncol(X);
  if (length(colnames(X)) == p) {
    name = colnames(X)
  } else {
    name = paste("X", 1:p, sep = "")
  }
  types_code = rep(NA, p); types_code[types == "con"] = 0; types_code[types == "bin"] = 1
  types_code[types == "tru"] = 2; types_code[types == "ter"] = 3
  cp = combn(p, 2); cp.col = ncol(cp); n0 = n * (n - 1) / 2
  n_X = apply(X, 2, function(x) {n_x(x = x, n)})
  n_X_sqd = sqrt(n0 - n_X)
  btoa_cp = matrix(n_X_sqd[cp], nrow = 2)
  btoa2 = btoa_cp[1, ] * btoa_cp[2, ] / n0
  K_b = pcaPP::cor.fk(X)
  K_b.lower = K_b[lower.tri(K_b)]
  K_a.lower = K_b.lower * btoa2
  zratios = vector(mode = "list", length = p)
  types_cp = matrix(types_code[cp], nrow = 2)
  X_12 = as.matrix(X[ , types_code == 1 | types_code == 2]); zratios[types_code == 1 | types_code == 2] = colMeans(X_12 == 0)
  X_3 = as.matrix(X[ , types_code == 3]); zratios[types_code == 3] = mattolist(rbind(colMeans(X_3 == 0), 1 - colMeans(X_3 == 2)))
  zratios_cp = matrix(zratios[cp], nrow = 2)
  types_mirror = types_cp[1, ] < types_cp[2, ]
  types_cp[ , types_mirror] = rbind(types_cp[2, types_mirror], types_cp[1, types_mirror])
  zratios_cp[ , types_mirror] = rbind(zratios_cp[2, types_mirror], zratios_cp[1, types_mirror])
  combs_cp = paste0(types_cp[1, ], types_cp[2, ]); combs = unique(combs_cp); R.lower = rep(NA, cp.col)
  for (comb in combs) {
    comb_select = combs_cp == comb; comb_select.len = sum(comb_select)
    bound_comb = get(paste("bound", comb, sep = "_"))
    cutoff = rep(NA, comb_select.len); K = K_a.lower[comb_select]
    zratio1 = matrix(unlist(zratios_cp[1, comb_select]), ncol = comb_select.len)
    zratio2 = matrix(unlist(zratios_cp[2, comb_select]), ncol = comb_select.len)
    if (method == "original") {
      R.lower[comb_select] = r_sol(K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = tol)
    } else if (method == "ml") {
      R.lower[comb_select] = r_ml(K = K, zratio1 = zratio1, zratio2 = zratio2, bound_comb = bound_comb, comb = comb)
    } else {
      cutoff = abs(K) > ratio * bound_comb(zratio1 = zratio1, zratio2 = zratio2)
      R.lower[comb_select][cutoff] = r_sol(K = K[cutoff], zratio1 = zratio1[ , cutoff], zratio2 = zratio2[ , cutoff], comb = comb, tol = tol)
      R.lower[comb_select][!(cutoff)] = r_ml(K = K[!(cutoff)], zratio1 = zratio1[ , !(cutoff)], zratio2 = zratio2[ , !(cutoff)], bound_comb = bound_comb, comb = comb)
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
  if (plotR) {
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

mattolist = function(X) {lapply(seq(ncol(X)), function(i) X[,i])}

r_sol = function(K, zratio1, zratio2, comb, tol) {
  bridge_comb = get(paste("bridge", comb, sep = "_")); K.len = length(K); out = rep(NA, K.len)
  zratio1 = as.matrix(zratio1); zratio2 = as.matrix(zratio2)
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

r_ml = function(K, zratio1, zratio2, bound_comb, comb) {
  ipol_comb = get(paste("ipol", comb, sep = "_"))
  K.len = length(K); out = rep(NA, K.len)
  K = K / bound_comb(zratio1 = zratios[1, type_comb], zratio2 = zratios[2, type_comb])
  if (ncol(zratio1) == 2) {zratio1[1, ] = zratio1[1, ] / zratio1[2, ]}
  if (ncol(zratio2) == 2) {zratio2[1, ] = zratio2[1, ] / zratio2[2, ]}
  out = ipol_comb(c(k, zratio1, zratio2)) / 10^7
  return(out)
}
