#' Evaluate an interpolant in a point
#'
#' An interpolant is a function returned by \code{\link{ipol}} which has prespecified values in some
#' points, and which fills in between with some reasonable values.
#'
#' @name interpolant
#' @param x The argument of the function. A function of more then one variable takes a
#' vector. \code{x} can also be a matrix of column vectors.
#' @param threads The number of threads to use for evaluation. All  interpolants created by
#' \code{ipol} are parallelized. If given a matrix argument \code{x}, the vectors can
#' be evaluated in parallel.
#' @param ... Other parameters. Currently used for simplex linear interpolants with the logical argument
#' The \code{"multilinear"} interpolant also has the argument \code{blend=c("linear","cubic","sigmoid")} where a
#' blending function can be chosen.
#' @return A numeric. If more than one point was evaluated, a vector.
#' @author Simen Gaure
NULL
