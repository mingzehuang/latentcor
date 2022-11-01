rm(list=ls())
library(boot)
library(latentcor)
# Create objective statistics
statistic <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  latent_R <- latentcor(X=d)$R[lower.tri(diag(ncol(d)))]
  return(latent_R)
}
inputdata=gen_data(types=c("con","bin","ter","tru"))$X
# bootstrapping with 1000 replications (R=1000)
results <- boot(data=inputdata, statistic=statistic, R=1000)

# view results
results

dcol=ncol(inputdata)
CI_lower=CI_upper=matrix(NA,dcol,dcol)
len_lowtri=sum(lower.tri(CI_lower))
boot_lower=boot_upper=rep(NA, len_lowtri)
# get 95% confidence interval (alpha=0.05)
alpha=0.05
for (i in 1:len_lowtri){
  boot_CI=boot.ci(results, conf=1-alpha, type="bca", index=i)
  boot_lower[i]=boot_CI$bca[4]
  boot_upper[i]=boot_CI$bca[5]
}
#CI_lower is the lower bound of confidence interval for lower triangular of R
#CI_upper is the upper bound of confidence interval for lower triangular of R
CI_lower[lower.tri(CI_lower)]=boot_lower
CI_upper[lower.tri(CI_upper)]=boot_upper

