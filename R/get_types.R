get_types = function(X) {
  X = as.matrix(X);
  n = nrow(X); p = ncol(X); X_num = matrix(NA, n, p); types = rep(NA, p)
  for (i in 1:p) {
    X_fac[ , i] = factor(X[ , i])
    level = levels(X_fac[ , i])
    X_num[ , i] = as.numeric(X_fac[ , i]) - 1
    if (length(level) <= 1) {
      stop("No variation in ", i, "th variable (", i, "th column of input data).")
    } else if (length(level) == 2) {
      types[i] = "bin"
    } else if (length(level) == 3) {
      types[i] = "ter"
    } else {
      if (length(X_num[ , i] == min(X_num[ , i])) / length(X_num[ , i]) > 0.1) {
        types[i] = "tru"
      } else {
        types[i] = "con"
      }
    }
  }
}
