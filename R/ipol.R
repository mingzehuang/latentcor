#' Methods for creating multivariate interpolations on hypercubes (originally in chebpol R package implemented by Simen Gaure).
#'
#' The primary method is \code{\link{ipol}} which
#' dispatches to some other method.  All the generated
#' \link{interpolant}s accept as an argument a matrix of column
#' vectors. The generated functions also accept an argument
#' \code{threads=getOption('ipol.threads')} to utilize more than
#' one CPU if a matrix of column vectors is evaluated.  The option
#' \code{ipol.threads} is initialized from the environment variable
#' \code{IPOL_THREADS} upon loading of the package. It defaults to \code{1}.
#'
#' The interpolants are ordinary R-objects and can be saved with \code{save()} and loaded
#' later with \code{load()} or serialized/unserialized with other tools, just like any R-object.
#' However, they contain calls to functions in the package, and while the author will make efforts
#' to ensure that generated interpolants are compatible with future versions of \pkg{ipol},
#' I can issue no such absolute guarantee.
#'
#'
#' @name interpolation
#' @aliases interpolation
#' @seealso \link{ipol}, \link{interpolant}
#' @useDynLib latentcor, .registration=TRUE, .fixes='C_'
#'
NULL
#'
#' Create an interpolating function from given values. Several interpolation methods are
#' supported.
#'
#' @param val array or function. Function values on a grid.
#' @param grid list. Each element is a vector of ordered grid-points for a
#' dimension.
#' @param ... Further arguments to the function, if \code{is.function(val)}. And some
#' extra arguments for interpolant creation described in section Details.
#' @return A \code{function(x, threads=getOption('chebpol.threads'))} defined on a hypercube, an \link{interpolant}
#' for the given function. The argument \code{x} can be a matrix of column
#' vectors which are evaluated in parallel in a number of threads.  The
#' function yields values for arguments outside the hypercube as well, though
#' it will typically be a poor approximation.  \code{threads} is an integer
#' specifying the number of parallel threads which should be used when
#' evaluating a matrix of column vectors.
#' @import stats geometry
#' @author Simen Gaure
#' @examples
#'
#' ## evenly spaced grid-points
#' su <- seq(0,1,length.out=10)
#' ## irregularly spaced grid-points
#' s <- su^3
#' ## create approximation on the irregularly spaced grid
#' ml1 <- ipol(exp(s), grid=list(s))
#' ## test it, since exp is convex, the linear approximation lies above
#' ## the exp between the grid points
#' ml1(su) - exp(su)
#'
#' ## multi dimensional approximation
#' f <- function(x) 10/(1+25*mean(x^2))
#' # a 3-dimensional 10x10x10 grid, first and third coordinate are non-uniform
#' grid <- list(s, su, sort(1-s))
#'
#' # make multilinear spline.
#' ml2 <- ipol(array(apply(expand.grid(grid), 1, f), c(10, 10, 10)), grid=grid)
#' # make 7 points in R3 to test them on
#' m <- matrix(runif(3*7),3)
#' rbind(true=apply(m,2,f), ml=ml2(m))
#'
#' @export
ipol <- function(val, grid=NULL, ...) {
  args <- list(...)
           if(is.null(grid)) stop('grid must be specified for multi linear interpolation')
           if(!is.list(grid)) grid <- list(grid)
           grid <- lapply(grid,as.numeric)
           if(unsortedgrid(grid)) stop("grid must be distinct ordered values")
           blend <- args[['blend']]
           if(is.null(blend)) return(mlappx.real(val,grid,...))
           newarg <- args[-match(c('blend'),names(args), nomatch=0)]
           return(blenddef(do.call(mlappx.real,c(list(val,grid), newarg)),blend))
}

unsortedgrid <- function(g) {
  any(sapply(g,function(s) is.unsorted(s,strictly=TRUE) && is.unsorted(-s,strictly=TRUE)))
}

