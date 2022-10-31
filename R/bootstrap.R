# Bootstrap 95% CI
statistic <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  latentcorr <- latentcor(X=d)$R[2,1]
  return(latentcorr)
}
# bootstrapping with 1000 replications
results <- boot(data=gen_data()$X, statistic=statistic, R=1000)

# view results
results
plot(results)

# get 95% confidence interval
boot.ci(results, type="bca")
