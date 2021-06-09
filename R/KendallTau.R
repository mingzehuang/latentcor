#' Kendall's tau correlation
#'
#' Calculate Kendall's tau correlation.
#' \deqn{ \hat{\tau}_{jk} = \frac{2}{n(n-1)}\sum_{1\le i<i'\le n} sign(X_{ji}-X_{ji'}) sign(X_{ki}-X_{ki'}) }
#' @param X A numeric matrix (n by p1).
#' @param Y A numeric matrix (n by p2).
#'
#' @rdname KendallTau
#'
#' @return \code{Kendall_matrix(X)} returns a p1 by p1 matrix of Kendall's tau correlation coefficients. \code{Kendall_matrix(X, Y)} returns a p1 by p2 matrix of Kendall's tau correlation coefficients.
#'
#' @export
#' @importFrom pcaPP cor.fk
Kendall_matrix <- function(X, Y = NULL){ # X and Y are matrix.
  X <- as.matrix(X) # In case that X is vector, this line gives you one column matrix.
  p1 <- ncol(X)
  if (is.null(Y)) {
    tau <- matrix(1, p1, p1)
    if (p1 > 1){
      for (i in 1:(p1 - 1)){
        for (j in (i + 1):p1){ # calculate KendallTau elementwise.
          tau[i, j] <- tau[j, i] <- KendallTau(X[, i], X[, j])
        }
      }
    }
  } else {
    Y <- as.matrix(Y) # In case that X and Y are vectors, this line gives you one column matrix.
    p2 <- ncol(Y)
    if (p1 == 1 & p2 == 1){
      tau <- KendallTau(X, Y)
    } else {
      tau <- matrix(0, p1, p2)
      for (i in 1:p1){
        for (j in 1:p2){ # calculate KendallTau elementwise.
          tau[i, j] <- KendallTau(X[, i], Y[, j])
        }
      }
    }
  }
  return(tau)
}
