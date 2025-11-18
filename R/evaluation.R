#' @title Numerical evaluation for different estimation methods.
#' @description Speed and accuracy comparison of two different estimation methods.
#' @param genfun A data generation function.
#' @param estfun_1 A function for first estimation method.
#' @param estfun_2 A function for second estimation method.
#' @param grid_list A list for grid points to be evaluated (each element of list is a vector represents ticklabels on a dimension). The number of list elements are the dimension of function inputs.
#' @param nrep Number of replications in simulation.
#' @param cores The numbers of cores (threads) of your machine to conduct parallel computing.
#' @param showplot Logical indicator. \code{showplot = TRUE} generates the heatmaps of output arrays. NULL if \code{showplot = FALSE}.
#' @param ... Other inputs for data generation or estimation functions to be passed through.
#' @return \code{evaluation} returns
#' \itemize{
#'      \item meanAE_1: An array for mean absolute error of first estimation method.
#'      \item meanAE_2: An array for mean absolute error of second estimation method.
#'      \item medianAE_1: An array for median absolute error of first estimation method.
#'      \item medianAE_2: An array for median absolute error of second estimation method.
#'      \item maxAE_1: An array for maximum absolute error of first estimation method.
#'      \item maxAE_2: An array for maximum absolute error of second estimation method.
#'      \item meanAE_diff: An array for mean absolute error of difference between two estimations.
#'      \item medianAE_diff: An array for median absolute error of difference between two estimations.
#'      \item maxAE_diff: An array for maximum absolute error of difference between two estimations.
#'      \item mediantime_1: An array for median time of first estimation method.
#'      \item mediantime_2: An array for median time of second estimation method.
#'      \item plot_meanAE_1: A plot for mean absolute error of first estimation method.
#'      \item plot_meanAE_2: A plot for mean absolute error of second estimation method.
#'      \item plot_medianAE_1: A plot for median absolute error of first estimation method.
#'      \item plot_medianAE_2: A plot for median absolute error of second estimation method.
#'      \item plot_maxAE_1: A plot for maximum absolute error of first estimation method.
#'      \item plot_maxAE_2: A plot for maximum absolute error of second estimation method.
#'      \item plot_meanAE_diff: A plot for mean absolute error of difference between two estimations.
#'      \item plot_medianAE_diff: A plot for median absolute error of difference between two estimations.
#'      \item plot_maxAE_diff: A plot for maximum absolute error of difference between two estimations.
#'      \item plot_mediantime_1: A plot for median time of first estimation method.
#'      \item plot_mediantime_2: A plot for median time of second estimation method.
#' }
#' @import doRNG
#' @importFrom stats median
#' @importFrom microbenchmark microbenchmark
#' @export

evaluation = function(genfun, estfun_1, estfun_2, grid_list, nrep = 100, showplot = FALSE, cores = detectCores(), ...) {
  grid_all = expand.grid(grid_list)
  registerDoFuture()
  plan(multicore, workers = cores)
  j = NULL
  value_vector =
    foreach (j = 1:nrow(grid_all), .combine = rbind) %dorng% {
      grid_input = as.numeric(grid_all[j, ])
      time_reps_1 = time_reps_2 = estimate_reps_1 = estimate_reps_2 = diff_reps = rep(NA, nrep)
      for (r in 1:nrep) {
        estimate_1 = estimate_2 = NULL
        simdata = genfun(grid_input, ...)
        time_reps_1[r] = median(microbenchmark::microbenchmark(estimate_1 <- estfun_1(simdata, ...), times = 5)$time)
        time_reps_2[r] = median(microbenchmark::microbenchmark(estimate_2 <- estfun_2(simdata, ...), times = 5)$time)
        estimate_reps_1[r] = estimate_1; estimate_reps_2[r] = estimate_2
        diff_reps[r] = estimate_1 - estimate_2
      }
      mean_AE_1 = mean(abs(estimate_reps_1 - grid_input[1])); mean_AE_2 = mean(abs(estimate_reps_2 - grid_input[1]))
      median_AE_1 = median(abs(estimate_reps_1 - grid_input[1])); median_AE_2 = median(abs(estimate_reps_2 - grid_input[1]))
      max_AE_1 = max(abs(estimate_reps_1 - grid_input[1])); max_AE_2 = max(abs(estimate_reps_2 - grid_input[1]))
      mean_AE_diff = mean(abs(diff_reps)); median_AE_diff = median(abs(diff_reps)); max_AE_diff = max(abs(diff_reps))
      median_time_1 = median(time_reps_1); median_time_2 = mean(time_reps_2)
      value_list = c(mean_AE_1, mean_AE_2, median_AE_1, median_AE_2, max_AE_1, max_AE_2, mean_AE_diff, median_AE_diff, max_AE_diff, median_time_1, median_time_2)
    }
  d_grid = length(grid_list)
  dim_value = NULL
  for (i in 1:d_grid) {
    dim_value = c(dim_value, length(grid_list[[i]]))
  }
  meanAE_1 = array(value_vector[ , 1], dim = dim_value)
  meanAE_2 = array(value_vector[ , 2], dim = dim_value)
  medianAE_1 = array(value_vector[ , 3], dim = dim_value)
  medianAE_2 = array(value_vector[ , 4], dim = dim_value)
  maxAE_1 = array(value_vector[ , 5], dim = dim_value)
  maxAE_2 = array(value_vector[ , 6], dim = dim_value)
  meanAE_diff = array(value_vector[ , 7], dim = dim_value)
  medianAE_diff = array(value_vector[ , 8], dim = dim_value)
  maxAE_diff = array(value_vector[ , 9], dim = dim_value)
  mediantime_1 = array(value_vector[ , 10], dim = dim_value)
  mediantime_2 = array(value_vector[ , 11], dim = dim_value)
  plot_meanAE_1 = plot_meanAE_2 = plot_medianAE_1 = plot_medianAE_2 = plot_maxAE_1 = plot_maxAE_2 =
  plot_meanAE_diff = plot_medianAE_diff = plot_maxAE_diff = plot_mediantime_1 = plot_mediantime_2 = NULL
  if (showplot) {
    for (i in c("meanAE_1", "meanAE_2", "medianAE_1", "medianAE_2", "maxAE_1", "maxAE_2", "meanAE_diff", "medianAE_diff", "maxAE_diff", "mediantime_1", "mediantime_2")) {
      eval_data = data.frame(get(i))
      rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
      assign(paste("plot", i, sep = "_"), heatmaply(eval_data, dendrogram = "none", main = i, margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
                                    grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), paste0(i, ":"))))
    }
  }
  return (list(meanAE_1 = meanAE_1, meanAE_2 = meanAE_2, medianAE_1 = medianAE_1, medianAE_2 = medianAE_2,
               maxAE_1 = maxAE_1, maxAE_2 = maxAE_2, meanAE_diff = meanAE_diff, medianAE_diff = medianAE_diff, maxAE_diff = maxAE_diff,
               mediantime_1 = mediantime_1, mediantime_2 = mediantime_2,
               plot_meanAE_1 = plot_meanAE_1, plot_meanAE_2 = plot_meanAE_2, plot_medianAE_1 = plot_medianAE_1, plot_medianAE_2 = plot_medianAE_2,
               plot_maxAE_1 = plot_maxAE_1, plot_maxAE_2 = plot_maxAE_2, plot_meanAE_diff = plot_meanAE_diff, plot_medianAE_diff = plot_medianAE_diff,
               plot_maxAE_diff = plot_maxAE_diff, plot_mediantime_1 = plot_mediantime_1, plot_mediantime_2 = plot_mediantime_2))
}
