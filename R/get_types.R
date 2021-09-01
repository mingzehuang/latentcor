#' Automatically determine types of each variable (continuous/binary/ternary/truncated) in a data matrix.
#'
#' @param X A numeric data matrix (n by p), where n is number of samples, and p is number of variables. Missing values (NA) are allowed.
#' @return \code{get_types} returns
#' \itemize{
#'       \item{types: }{A vector of length p indicating the type of each of the p variables in \code{X}. Each element must be one of \code{"con"} (continuous), \code{"bin"} (binary), \code{"ter"} (ternary) or \code{"tru"} (truncated).}
#' }
#' @export
#' @examples
#' X = gen_data(types = c("ter", "con"))$X
#' get_types(X = X)

get_types = function(X) {
  X = data.matrix(X)
  if (!(is.numeric(X))) {
    stop("Input data matrix should be numeric.")
  }
  # Get the number of variables
  p = ncol(X); types = rep(NA, p)
  # Determine the type of each variable
  for (i in 1:p) {
    # Extract the number of unique levels
    level = unique(X[ , i])
    if (length(level) <= 1) {
      # Only one level means no variation
      stop("No variation in ", i, "th variable (", i, "th column of input data).")
    } else if (length(level) == 2) {
      # Two levels means binary
      types[i] = "bin"
    } else if (length(level) == 3) {
      # Three levels means ternary
      types[i] = "ter"
    } else if (length(level) > 3 & length(level) <= 10) {
      # Since ordinals of more than 3 levels are not supported, alert the user those are treated as continuous
      message("ordinal levels between 4 and 10 will be approximated by continuous type.")
      types[i] = "con"
    } else {
      if (mean(X[ , i] == min(X[ , i], rm.na = TRUE)) > 0.15) {
        types[i] = "tru"
      } else {
        types[i] = "con"
      }
    }
  }
  return(types)
}
