#' Mixed type simulation data generator
#' \code{GenData} is used to generate two sets of data of mixed types for sparse CCA under the Gaussian copula model.
#' @param n Sample size
#' @param copula1 Copula type for the first dataset. U1 = f(Z1), which could be either "exp", "cube".
#' @param copula2 Copula type for the second dataset. U2 = f(Z2), which could be either "exp", "cube".
#' @param type1 Type of the first dataset \code{X1}. Could be "continuous", "trunc", "binary", "ternary".
#' @param type2 Type of the second dataset \code{X2}. Could be "continuous", "trunc", "binary", "ternary".
#' @param rho True correlation for auto correlated data.
#' @param muZ Mean of latent multivariate normal.
#' @param Sigma True correlation matrix of latent variable \code{Z1} (p1 by p1).
#' @param p1 Dimension of variables belong to \code{type1}.
#' @param p2 Dimension of variables belong to \code{type2}.
#' @param c1 Constant threshold for \code{X1} needed for "trunc", "binary", "ternary" and "ordinal" data type - the default is NULL.
#' @param c2 Constant threshold for \code{X2} needed for "trunc", "binary", "ternary" and "ordinal" data type - the default is NULL.
#' @param q1 quantile threshold for \code{X1} needed for "trunc", "binary", "ternary" and "ordinal" data type - the default is NULL.
#' @param q2 quantile threshold for \code{X2} needed for "trunc", "binary", "ternary" and "ordinal" data type - the default is NULL.
#' @return \code{GenData} returns a list containing
#' \itemize{
#'       \item{Z1: }{latent numeric data matrix (n by p1).}
#'       \item{Z2: }{latent numeric data matrix (n by p2).}
#'       \item{X1: }{observed numeric data matrix (n by p1).}
#'       \item{X2: }{observed numeric data matrix (n by p2).}
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
#'
GenData = function(n, type1, type2, p1 = 1, p2 = 1, rho = NULL, copula1 = NULL, copula2 = NULL, muZ = NULL, Sigma = NULL, c1 = NULL, c2 = NULL, q1 = NULL, q2 = NULL){
  if (is.null(muZ)) {muZ = rep(0, p1 + p2)}
  if (is.null(Sigma) & is.null(rho)) {
    stop("correlation rho or covariance matrix Sigma needs to be specified.")
  } else if (is.null(Sigma)) {
    Sigma = autocor(p1 + p2, rho)
    }
  dat = MASS::mvrnorm(n, mu = muZ, Sigma = Sigma) # generate a data matrix of size.
  Z1 = as.matrix(dat[ , 1:p1]); Z2 = as.matrix(dat[ , (p1 + 1):(p1 + p2)])
  X1 = fromZtoX(Z = Z1, copula = copula1, type = type1, q = q1, c = c1)
  X2 = fromZtoX(Z = Z2, copula = copula2, type = type2, q = q2, c = c2)
  return(list(muZ = muZ, Sigma = Sigma, Z1 = Z1, Z2 = Z2, X1 = X1, X2 = X2))
}


#' Plot latent correlation
#' \code{LatentPlot} is to plot heatmap for latent correlation matrix.
#' @param latentcor latent correlation matrix.
#' @param xlab label for horizontal axis.
#' @param ylab label for vertical axis.
#' @param main main title for graphs.
#' @import heatmaply hrbrthemes ggplot2
#' @importFrom utils data
#' @return \code{LatentPlot} returns a heatmap plot for latent correlation matrix.
#' @example man/examples/estimateR_ex.R
#' @export
LatentPlot = function(latentcor, xlab = "Variable 1", ylab = "Variable 2", main = "Latent Correlation for Variable 1 vs. Variable 2") {
  out = heatmaply(latentcor, dendrogram = "none",
                       xlab = xlab, ylab = ylab,
                       main = main,
                       margins = c(60,100,40,20),
                       grid_color = "white",
                       grid_width = 0.00001,
                       titleX = TRUE,
                       hide_colorbar = FALSE,
                       branches_lwd = 0.1,
                       label_names = c(paste0(ylab, ":"), paste0(xlab, ":"), "Latent Correlation:"),
                       fontsize_row = 10, fontsize_col = 10,
                       labCol = colnames(data),
                       labRow = rownames(data),
                       heatmap_layers = theme(axis.line=element_blank()))
  return(out)
}
