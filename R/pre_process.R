#' Pre-process of raw data to deal with non-numeric variable, then get suggested variable types.
#'
#' @param raw_data A data frame or data matrix (n by p), where n is number of samples, and p is number of variables. Missing values (NA) are allowed, in which case the estimation is based on pairwise complete observations. Non-numeric variables are allowed, in which case users should provide ordinal orders of non-numeric values by \code{ordinals}.
#' @param ordinals A list to convert characters into ordinal numbers. Each element of the list is corresponding to each variable. For the numeric variable, please put NA; for the non-numeric variable, please put the character in ascending order to be converted to ordinal variable. The default value is NULL for the case that all variables are numeric. See vignette for detail.
#' @return \code{pre_process} returns
#' \itemize{
#'       \item{X: }{A numeric data matrix (n by p), where n is number of samples, and p is number of variables. Missing values (NA) are allowed, in which case the estimation is based on pairwise complete observations.}
#'       \item{types: }{A vector of length p indicating the type of each of the p variables in \code{X}. Each element must be one of \code{"con"} (continuous), \code{"bin"} (binary), \code{"ter"} (ternary), \code{"tru"} (truncated) or \code{"utr"} (upper truncated). If the vector has length 1, then all p variables are assumed to be of the same type that is supplied. The default value is \code{"con"} which means all variables are continuous.}
#' }
#' @export
#' @examples
#' raw_data = data.frame(matrix(NA, 3, 3))
#' raw_data$X1 = 1:3
#' raw_data$X2 = c("medium", "small", "large")
#' raw_data$X3 = 7:9
#' pre_process = pre_process(raw_data = raw_data, ordinals = list(NA, c("small", "medium", "large")))
#' pre_process$X
#' pre_process$types
#'

pre_process = function(raw_data, ordinals = NULL) {
  X = data.frame(raw_data); n = nrow(X); p = ncol(X); types = rep(NA, p)
  for (i in 1:p) {
    if (!(is.numeric(X[ , i]))) {
      X[ , i] = ord(X[ , i], ordinals[[i]])
    }
    level = unique(X[ , i])
    if (length(level) <= 1) {
      stop("No variation in ", i, "th variable (", i, "th column of input data).")
    } else if (length(level) == 2) {
      types[i] = "bin"
    } else if (length(level) == 3) {
      types[i] = "ter"
    } else if (length(level) > 3 & length(level) <= 10) {
      message("ordinal levels between 4 and 10 will be approximated by continuous type.")
      types[i] = "con"
    } else {
      if (length(X[ , i] == min(X[ , i], rm.na = TRUE)) / length(X[ , i]) > 0.15) {
        types[i] = "tru"
      } else if (length(X[ , i] == max(X[ , i], rm.na = TRUE)) / length(X[ , i]) > 0.15) {
        types[i] = "utr"
      } else {
        types[i] = "con"
      }
    }
  }
  return(list(X = as.matrix(X), types = types))
}
