#' Multilinear interpolation on a grid
#'
#' Multilinear interpolation on an arbitrary Cartesian product.
#'
#' A call \code{fun <- mlappx(val,grid)} creates a multilinear interpolant on
#' the grid.  The value on the grid points will be exact, the value between the
#' grid points is a convex combination of the values in the corners of the
#' hypercube surrounding it.
#'
#' If \code{val} is a function it will be evaluated on the grid.
#'
#' @aliases mlappx mlappxf
#' @param val Array or function. Function values on a grid, or the function
#' itself. If it is the values, the \code{dim}-attribute must be appropriately
#' set.
#' @param grid A list.  Each element is a vector of ordered grid-points for a
#' dimension.  These need not be Chebyshev-knots, nor evenly spaced.
#' @param ... Further arguments to the function, if \code{is.function(val)}.
#' @return A \code{function(x)} defined on the hypercube, approximating the
#' given function.  The function yields values for arguments outside the
#' hypercube as well, as a linear extension.
#' @keywords internal
mlappx.real <- function(val, grid, ...) {
  if(is.numeric(grid)) grid <- list(grid)
  if(any(sapply(grid,is.unsorted))) {
    if(!is.function(val)) stop('Grid points must be ordered in increasing order')
    grid <- lapply(grid,sort)
  }
  gl <- prod(sapply(grid,length))
  if(length(val) != gl)
    stop("length of values ",length(val)," do not match size of grid ",gl)
  val <- as.numeric(val)
  vectorfun(function(x,threads=getOption('latentcor.threads'),
                     blend=c('linear','cubic','sigmoid','parodic','square')) {
    blend <- switch(match.arg(blend),linear=0L,sigmoid=1L,parodic=2L,cubic=3L,square=4L)
    .Call(C_evalmlip,grid,val,x,threads,blend)
  },
  arity=length(grid),
  domain=lapply(grid,range))
}
