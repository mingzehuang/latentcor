#' Plot true correlation vs estimated correlation
#' \code{PlotCompare} is to check unbiasness of estimation by plotting true correlation from simulation data vs estimated correlation from simulation data.
#' @param pairlist list for data pairs.
#' @param namelist list for names of data pairs.
#' @import ggplot2
#' @return \code{PlotCompare} returns a plot for data1 against data2 and 45 degree benchmark line.
#' @example man/examples/estimateR_ex.R
#' @export
PlotCompare <- function(pairlist, namelist) {
  l_list <- length(pairlist)
  for (i in 1:l_list) {
    df <- data.frame(pairlist[[i]])
    colnames(df) = namelist[[i]]
    ggplot(df, aes(df[paste(namelist[[i]][1])], df[paste(namelist[[i]][2])]), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
  }
}
