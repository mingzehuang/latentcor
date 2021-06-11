
fromZtoX = function(Z, copula, type, c) {
  if(copula == "exp"){
    Z = exp(Z)
  }else if(copula == "cube"){
    Z = Z^3
  }
  if(type == "continuous") {
    X = Z
  } else if (ncol(c) != ncol(Z)) {
    stop("The length of threshold vector c does not match with the dimension of the data p.")
  } else if(type == "trunc") {
    X = ifelse(Z > matrix(c, nrow = nrow(Z), ncol = ncol(Z), byrow = T), Z, 0)
  } else if (type == "binary") {
    X = ifelse(Z > matrix(c, nrow = nrow(Z), ncol = ncol(Z), byrow = T), 1, 0)
  } else if (type == "ternary") {
    X = ifelse(Z >= matrix(apply(c, 2, max), nrow = nrow(Z), ncol = ncol(Z), byrow = T), 2, 1)
    X[Z <= matrix(apply(c, 2, min), nrow = nrow(Z), ncol = ncol(Z), byrow = T)] = 0
  }
  return(X)
}
