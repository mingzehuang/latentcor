#' Plot true correlation vs estimated correlation
#' \code{PlotCompare} is to check unbiasness of estimation by plotting true correlation from simulation data vs estimated correlation from simulation data.
#' @param pairlist list for data pairs.
#' @param namelist list for names of data pairs.
#' @import ggplot2
#' @return \code{PlotCompare} returns a plot for data1 against data2 and 45 degree benchmark line.
#' @example man/examples/estimateR_ex.R
#' @export
PlotCompare <- function(pairlist, namelist, title) {
  l_list <- length(pairlist)
  plotlist <- list(rep(NA, l_list))
  for (i in 1:l_list) {
    df <- data.frame(pairlist[[i]])
    colnames(df) = namelist[[i]]
    print(ggplot(df, aes(x = get(paste(namelist[[i]][1])), y = get(paste(namelist[[i]][2]))))
          + geom_point(color = "blue") + geom_abline(intercept = 0, slope = 1, color = "red")
          +ggtitle(title) +
            xlab(paste(namelist[[i]][1])) + ylab(paste(namelist[[i]][2])))
  }
}
