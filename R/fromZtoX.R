autocor = function(p, rho){
  if (abs(rho) > 0.99){stop("correlation rho must be between -0.99 and 0.99.")}
  return(rho^abs(outer(1:p, 1:p, "-")))
}

fromZtoX = function(Z, copula, type, c, q) {
  p = ncol(Z)
  if(copula == "exp"){
    Z = exp(Z)
  }else if(copula == "cube"){
    Z = Z^3
  }
  if(type == "continuous") {
    X = Z
  } else {
    if (is.null(q) & is.null(c)) {
      stop("quantile threshold q or level threshold c needs to be defined for truncated, binary, ternary data type.")
    } else if (is.null(c)) {
      q = as.matrix(q)
      q = matrix(q, nrow = nrow(q), ncol = p)
      for (i in 1:p) {c = cbind(quantile(Z[ , i], q[ , i]))}
    } else {
      c = as.matrix(c)
      c = matrix(c, nrow = nrow(c), ncol = p)
    }
    if (type == "binary") {
    X = ifelse(Z > matrix(c, nrow = nrow(Z), ncol = ncol(Z), byrow = T), 1, 0)
    } else if(type == "trunc") {
      X = ifelse(Z > matrix(c, nrow = nrow(Z), ncol = ncol(Z), byrow = T), Z, 0)
    } else if (type == "ternary") {
      X = ifelse(Z >= matrix(apply(c, 2, max), nrow = nrow(Z), ncol = ncol(Z), byrow = T), 2, 1)
      X[Z <= matrix(apply(c, 2, min), nrow = nrow(Z), ncol = ncol(Z), byrow = T)] = 0
    }
  }
  return(X)
}
