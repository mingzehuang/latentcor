#' Plot true correlation vs estimated correlation
#' \code{PlotCompare} is to check unbiasness of estimation by plotting true correlation from simulation data vs estimated correlation from simulation data.
#' @param data1 First data series.
#' @param data2 Second data series.
#' @param name1 Name for first data series.
#' @param name2 Name for second data series.
#' @import ggplot2
#' @return \code{PlotCompare} returns a plot for data1 against data2 and 45 degree benchmark line.
#' @example man/examples/estimateR_ex.R
#' @export
PlotCompare <- function(data1, data2, name1, name2) {
df <- data.frame(c(data1), c(data2))
colnames(df) = c(name1, name2)
  out <- ggplot(df, aes(df[paste(name1)], df[paste(name2)]), color = "blue") + geom_point() + geom_abline(intercept = 0, slope = 1, color="red")
return(out)
}
