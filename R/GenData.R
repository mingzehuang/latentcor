#'
#' @title Mixed type simulation data generator
#' @description \code{GenData} is used to generate data of mixed types for the Gaussian copula model. It can generate single, pair or multiple variables.
#' @param n a positive interger to specify sample size. The default value is 100.
#' @param types A length of p vector as types for each variable. Could be "con" (continuous), "bin" (binary), "tru" (truncated) or "ter" (ternary).
#'              The number of variables p is determined based on supplied length of types (or something along those lines).
#'              The default value is c("ter", "con") which specifies a pair of variables, the first one is ternary, the second one is continuous.
#' @param rhos A vector with length of 1 as true correlation for between any two variables (e.g. off-diagonal elements of correlation matrix, so that the correlation matrix is \code{matrix(c(1, rhos, rhos, 1), 2, 2)}) or length of \code{sum(lower.tri(diag(p)))} as lower triangular elements of correlation matrix (e.g. \code{rhos = c(.3, .5, .7)} means the correlation matrix is \code{matrix(c(1, .3, .5, .3, 1, .7, .5, .7, 1), 3, 3)}). The default value is 0.5 which means correlations between any two variables are 0.5.
#' @param copulas A vector with length of 1 for copula type of all variables or p for copula types for each variable, e.g. U = f(Z). Which could be "no", "expo" (exponential) or "cube" (cube).
#'                The default value is "no" which means no copula transformation for any variables.
#' @param XP A length of p list of proportions of zeros (and ones for ternary) corresponding to each variable.
#'           If null, the following values are automatically generated as elements of \code{XP} list for the corresponding data types:
#'           For continuous variable, the corresponding proportion should be NA;
#'           for binary or truncated variable, the corresponding proportion should be a number between 0 and 1 represents the proportion of zeros;
#'           for ternary variable, the corresponding proportion should be a pair of numbers between 0 and 1, the first number indicates the the proportion of zeros, the second number indicates the proportion of ones. The sum of a pair of numbers should between 0 and 1.
#'           The default value for continuous variable is NA, for binary/truncated variable is 0.5, for ternary variable is c(0.3, 0.5).
#'           Otherwise, the list element should be NA if the corresponding variable is continuous; a number between 0 and 1 if the corresponding variable is binary/truncated; a vector of two numbers between 0 and 1 if the corresponding variable is ternary (first number smaller then the second, sum of them smaller than 1).
#'
#' @param showplot TRUE or FALSE if you want to plot data if numbers of variables (features) no more than 3. The default value is FALSE, not to generate plot.
#' @return \code{GenData} returns a list containing
#' \itemize{
#'       \item{X: }{observed numeric data matrix (n by p).}
#'       \item{plotX: }{Visualize data matrix X to check if it's consistent with types and latent correlations settings.
#'                     Histogram of data X if it only has one variable. 2D Scatter plot if it has 2 variables. 3D scatter plot if it has 3 variables.}
#' }
#' @export
#' @importFrom stats rnorm
#' @importFrom MASS mvrnorm
#' @importFrom graphics hist
#' @importFrom plotly plot_ly layout
#' @examples
#' # Generate single continuous variable with exponential transformation (always greater than 0)
#' # and show histogram.
#' simdata = GenData(n = 100, copulas = "expo", types = "con", showplot = FALSE)
#' X = simdata$X; plotX = simdata$plotX
#' # Generate a pair of variables (ternary and continuous) with default proportions
#' # and without copula transformation.
#' simdata = GenData()
#' X = simdata$X
#' # Generate 3 variables (binary, ternary and truncated)
#' # corresponding copulas for each variables are "no" (no transformation),
#' # "cube" (cube transformation) and "cube" (cube transformation).
#' # binary variable has 30% of zeros, ternary variable has 20% of zeros
#' # and 40% of ones, truncated variable has 50% of zeros.
#' # Then show the 3D scatter plot (data points project on either 0 or 1 on Axis X1;
#' # on 0, 1 or 2 on Axas X2; on positive domain on Axis X3)
#' simdata = GenData(n = 100, rhos = c(.3, .4, .5), copulas = c("no", "cube", "cube"),
#'           types = c("bin", "ter", "tru"), XP = list(.3, c(.2, .4), .5), showplot = TRUE)
#' X = simdata$X; plotX = simdata$plotX
#' # Check the proportion of zeros for the binary variable.
#' sum(simdata$X[ , 1] == 0)
#' # Check the proportion of zeros and ones for the ternary variable.
#' sum(simdata$X[ , 2] == 0); sum(simdata$X[ , 2] == 1)
#' # Check the proportion of zeros for the truncated variable.
#' sum(simdata$X[ , 3] == 0)

GenData = function(n = 100, types = c("ter", "con"), rhos = .5, copulas = "no", XP = NULL, showplot = FALSE) {
  if (length(n) != 1 | n <= 0) {stop("n should be a positive interger as sample size")}
# types = match.arg(types, c("con", "bin", "tru", "ter", "qua", "qui", "sen", "sep", "oct", "nov", "den", "dtr"), several.ok = TRUE)
  types = match.arg(types, c("con", "bin", "tru", "ter"), several.ok = TRUE)
  p = length(types)
  copulas = match.arg(copulas, c("no", "expo", "cube"), several.ok = TRUE); p.copulas = length(copulas)
  if (p.copulas == 1) {
    copulas = rep(copulas, p)
  } else if (p.copulas != p) {
    stop("copulas should have the same length as types, so that each copula corresponds to a variable (feature).")
  }
  if (is.null(XP)) {
    XP = vector(mode = "list", length = p)
    for (type in unique(types)) {
      XP[types == type] = switch(type, "con" = NA, "bin" = .5, "tru" = .5, "ter" = list(c(.3, .5))
                                 # , "qua" = list(c(.2, .2, .2)),
                                 # "qui" = list(c(.2, .2, .2, .2)), "sen" = list(c(.1, .1, .1, .1, .1)),
                                 # "sep" = list(c(.1, .1, .1, .1, .1, .1)), "oct" = list(c(.1, .1, .1, .1, .1, .1, .1)),
                                 # "nov" = list(c(.1, .1, .1, .1, .1, .1, .1, .1)), "den" = list(c(.1, .1, .1, .1, .1, .1, .1, .1, .1)),
                                 # "dtr" = list(c(.3, .5))
                                 )
    }
  } else if(!(is.list(XP)) | length(XP) != p) {
    stop("XP should be a list has the same length as types, so that each element is a set of proportion(s) corresponds to a variable (feature).")
  } else if (any(na.omit(unlist(XP)) <= 0) | any(na.omit(unlist(XP)) >= 1)) {
    stop("The proprotions should always between 0 and 1. Otherwise please consider to degenerate your data type.")
  }
  if (p == 1) {
    Z = as.matrix(rnorm(n = n))
  } else {
    if (length(copulas) == 1) {copulas = rep(copulas, p)}
    Sigma.lower = diag(0, p); lowerpart = lower.tri(Sigma.lower); rhos.len = length(rhos)
    if (rhos.len != 1 & rhos.len != sum(lowerpart)) {
      stop("Length of rhos should fit for lower triangular part of latent correplation matrix.")
    } else {
      Sigma.lower[lowerpart] = rhos; Sigma = Sigma.lower + t(Sigma.lower) + diag(1, p)
      Z = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    }
  }
  X = sapply(seq(p), function(i) {fromZtoX(z = Z[ , i], type = types[i], copula = copulas[i], xp = XP[[i]])})
  plotX = NULL
  if (p == 1 & showplot) {
    plotX = hist(X)
  } else if (p == 2 & showplot) {
    plotX = plot_ly(data.frame(X), x = ~X1, y = ~X2, type = "scatter", mode = "markers")
  } else if (p == 3 & showplot) {
    plotX = plot_ly(data.frame(X), x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", mode = "markers")
  }
  return(list(X = X, plotX = plotX))
}
