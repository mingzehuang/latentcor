

evaluation = function(genfun, estfun_1, estf_2, grid_list, nrep = 100, cores = detectCores(), ...) {
  grid_all = expand.grid(grid_list)
  registerDoFuture()
  plan(multicore, workers = cores)
  j = NULL
  value_vector =
    foreach (j = 1:nrow(grid_all), .combine = c) %dopar% {
      grid_input = as.numeric(grid_all[j, ])
      time_reps_1 = rep(NA, nrep); time_reps_2 = rep(NA, nrep); error_reps_1 = rep(NA, nrep); error_reps_2 = rep(NA, nrep); diff_reps = rep(NA, nrep)
      for (r in 1:nrep) {
        simdata = genfun(grid_input, ...)
        time_reps_1[r] = median(microbenchmark::microbenchmark(estimate_reps_1[r] = estfun_1(simdata, ...), times = 5)$time)
        time_reps_2[r] = median(microbenchmark::microbenchmark(estimate_reps_2[r] = estfun_2(simdata, ...), times = 5)$time)
        diff_reps[r] = estimate_reps_1[r] - estimate_reps_2[r]
      }
      mean_AE_1 = mean(abs(estimate_reps_1 - grid_input[1])); mean_AE_2 = mean(abs(estmate_reps_2 - grid_input[1]))
      median_AE_1 = median(abs(estimate_reps_1 - grid_input[1])); median_AE_2 = median(abs(estmate_reps_2 - grid_input[1]))
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
  mean_AE_1_array = array(value_vector[ , 1], dim = dim_value)
  mean_AE_2_array = array(value_vector[ , 2], dim = dim_value)
  median_AE_1_array = array(value_vector[ , 3], dim = dim_value)
  median_AE_2_array = array(value_vector[ , 4], dim = dim_value)
  max_AE_1_array = array(value_vector[ , 5], dim = dim_value)
  max_AE_2_array = array(value_vector[ , 6], dim = dim_value)
  mean_AE_diff = array(value_vector[ , 7], dim = dim_value)
  median_AE_diff = array(value_vector[ , 8], dim = dim_value)
  max_AE_diff = array(value_vector[ , 9], dim = dim_value)
  median_time_1 = array(value_vector[ , 10], dim = dim_value)
  median_time_2 = array(value_vector[ , 11], dim = dim_value)

  return (list(mean_AE_1_array = mean_AE_1_array, mean_AE_2_array = mean_AE_2_array, median_AE_1_array = median_AE_1_array, median_AE_2_array = median_AE_2_array,
               max_AE_1_array = max_AE_1_array, max_AE_2_array = max_AE_2_array, mean_AE_diff = mean_AE_diff, median_AE_diff = median_AE_diff, max_AE_diff = max_AE_diff,
               median_time_1 = median_time_1, median_time_2 = median_time_2))
}
