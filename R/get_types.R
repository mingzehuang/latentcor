#' Automatically determine types of each variable (continuous/binary/ternary/truncated) in a data matrix.
#'
#' @param X A numeric data matrix (n by p), where n is number of samples, and p is number of variables. Missing values (NA) are allowed.
#' @param tru_prop A scalar between 0 and 1 indicating the minimal proportion of zeros that should be present in a variable to be treated as \code{"tru"} (truncated type or zero-inflated) rather than as \code{"con"} (continuous type). The default value is 0.05 (any variable with more than 5\% of zero values among n samples is treated as truncated or zero-inflated)
#'
#'
#' @return \code{get_types} returns
#' \itemize{
#'       \item types: A vector of length p indicating the type of each of the p variables in \code{X}. Each element is one of \code{"con"} (continuous), \code{"bin"} (binary), \code{"ter"} (ternary) or \code{"tru"} (truncated).
#' }
#' @export
#' @examples
#' X = gen_data(types = c("ter", "con"))$X
#' get_types(X)
#'
get_types = function(X, tru_prop = 0.05) {
  X = as.matrix(X)
  if (!(is.numeric(X))) {
    stop("Input data matrix should be numeric.")
  }
  # Get the number of variables
  p = ncol(X); types = rep(NA, p)
  # Determine the type of each variable
  for (i in 1:p) {
    # Extract the number of unique levels
    x = X[ , i]
    levels = unique(x[!(is.na(x))])
    if (length(levels) <= 1) {
      # Only one level means no variation
      stop("No variation in ", i, "th variable (", i, "th column of input data).")
    } else if (length(levels) == 2) {
      # Two levels means binary
      types[i] = "bin"
    } else if (length(levels) == 3) {
      # Three levels means ternary
      types[i] = "ter"
    } else if (length(levels) > 3) {
      # More than 3 levels are detected - could be truncated or continuous
      if (length(levels) < 10){
        # Since ordinals of more than 3 levels are not supported, alert the user those are treated as continuous
        message("ordinal levels between 4 and 10 will be approximated by either continuous or truncated type.")
      }
      if ((min(x, na.rm = TRUE) == 0) & (mean(x == 0, na.rm = TRUE) > tru_prop)) {
        # If the minimal value is zero and there are at least tru_prop zeros -> truncated (zero-inflated)
        types[i] = "tru"
      } else {
        types[i] = "con"
      }
    }
  }
  return(types)
}
