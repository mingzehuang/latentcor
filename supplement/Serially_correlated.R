# All cases (AR(1))
rm(list = ls())
library(MASS)
library(mvtnorm)
library(ggplot2)
library(latentcor)
PlotPair = function(datapair, namepair = c("X", "Y"), title = "Plot X vs Y") {
  df = data.frame(datapair)
  colnames(df) = namepair
  print(ggplot(df, aes(x = datapair[ , 1], y = datapair[ , 2]))
        + geom_point(color = "blue") + geom_abline(intercept = 0, slope = 1, color = "red")
        +ggtitle(title) + xlab(namepair[1]) + ylab(namepair[2]) + xlim(-1, 1) + ylim(-1, 1))
}
corrrep = seq(-1, 1, by = .01)
Rrep = matrix(NA, length(corrrep), 10)
n = 500
rho = matrix(NA, n, n)
rho = abs(col(rho) - row(rho))
r = 0.9
rho = r^rho
p = 2
r_bar = (sum(rho)-n)/(n*(n-1))
for (i in 1:length(corrrep)) {
    Z = matrix(mvrnorm(mu = rep(0, n * p), Sigma = kronecker(rho, matrix(c(1, corrrep[i], corrrep[i], 1), p, p))), ncol = p, byrow = TRUE)
    latentcor_CC = latentcor(Z, types = c("con", "con"))$Rpointwise[2, 1]
    X = Z
    Z_1 = Z[ , 1]; Z_bar_1 = quantile(Z_1, .5); X[ , 1][Z_1 <= Z_bar_1] = 0; X[ , 1][Z_1 > Z_bar_1] = 1
    latentcor_BC = latentcor(X, types = c("bin", "con"))$Rpointwise[2, 1]
    Z_2 = Z[ , 2]; Z_bar_2 = quantile(Z_2, .5); X[ , 2][Z_2 <= Z_bar_2] = 0; X[ , 2][Z_2 > Z_bar_2] = 1
    latentcor_BB = latentcor(X, types = c("bin", "bin"))$Rpointwise[2, 1]
    X = Z
    X[ , 1][Z_1 <= Z_bar_1] = 0
    latentcor_TC = latentcor(X, types = c("tru", "con"))$Rpointwise[2, 1]
    X[ , 2][Z_2 <= Z_bar_2] = 0; X[ , 2][Z_2 > Z_bar_2] = 1
    latentcor_TB = latentcor(X, types = c("tru", "bin"))$Rpointwise[2, 1]
    X = Z
    X[ , 1][Z_1 <= Z_bar_1] = 0
    X[ , 2][Z_2 <= Z_bar_2] = 0
    latentcor_TT = latentcor(X, types = c("tru", "tru"))$Rpointwise[2, 1]
    X = Z
    Z_bar_11 = quantile(Z_1, .3); Z_bar_12 = quantile(Z_1, .8)
    X[ , 1] = 1; X[ , 1][Z_1 <= Z_bar_11] = 0; X[ , 1][Z_1 > Z_bar_12] = 2
    latentcor_NC = latentcor(X, types = c("ter", "con"))$Rpointwise[2, 1]
    X[ , 2][Z_2 <= Z_bar_2] = 0; X[ , 2][Z_2 > Z_bar_2] = 1
    latentcor_NB = latentcor(X, types = c("ter", "bin"))$Rpointwise[2, 1]
    X = Z
    X[ , 1] = 1; X[ , 1][Z_1 <= Z_bar_11] = 0; X[ , 1][Z_1 > Z_bar_12] = 2
    X[ , 2][Z_2 <= Z_bar_2] = 0
    latentcor_NT = latentcor(X, types = c("ter", "tru"))$Rpointwise[2, 1]
    Z_bar_21 = quantile(Z_2, .3); Z_bar_22 = quantile(Z_2, .8)
    X[ , 2] = 1; X[ , 2][Z_2 <= Z_bar_21] = 0; X[ , 2][Z_2 > Z_bar_22] = 2
    latentcor_NN = latentcor(X, types = c("ter", "ter"))$Rpointwise[2, 1]
    Rrep[i, 1] = latentcor_CC; Rrep[i, 2] = latentcor_BC; Rrep[i, 3] = latentcor_BB;
    Rrep[i, 4] = latentcor_TC; Rrep[i, 5] = latentcor_TB; Rrep[i, 6] = latentcor_TT;
    Rrep[i, 7] = latentcor_NC; Rrep[i, 8] = latentcor_NB; Rrep[i, 9] = latentcor_NT;
    Rrep[i, 10] = latentcor_NN
  }
plot_CC = PlotPair(datapair = cbind(corrrep, Rrep[ , 1]),
         namepair = c("True latent correlation", "Estimated latent correlation"),
         title = paste("continuous vs. continuous", " (autocorrelation = ", r, ")", sep = ""))
plot_BC = PlotPair(datapair = cbind(corrrep, Rrep[ , 2]),
         namepair = c("True latent correlation", "Estimated latent correlation"),
         title = paste("binary vs. continuous", " (autocorrelation = ", r, ")", sep = ""))
plot_BB = PlotPair(datapair = cbind(corrrep, Rrep[ , 3]),
         namepair = c("True latent correlation", "Estimated latent correlation"),
         title = paste("binary vs. binary", " (autocorrelation = ", r, ")", sep = ""))
plot_TC = PlotPair(datapair = cbind(corrrep, Rrep[ , 4]),
         namepair = c("True latent correlation", "Estimated latent correlation"),
         title = paste("truncated vs. continuous", " (autocorrelation = ", r, ")", sep = ""))
plot_TB = PlotPair(datapair = cbind(corrrep, Rrep[ , 5]),
         namepair = c("True latent correlation", "Estimated latent correlation"),
         title = paste("truncated vs. binary", " (autocorrelation = ", r, ")", sep = ""))
plot_TT = PlotPair(datapair = cbind(corrrep, Rrep[ , 6]),
         namepair = c("True latent correlation", "Estimated latent correlation"),
         title = paste("truncated vs. truncated", " (autocorrelation = ", r, ")", sep = ""))
plot_NC = PlotPair(datapair = cbind(corrrep, Rrep[ , 7]),
         namepair = c("True latent correlation", "Estimated latent correlation"),
         title = paste("ternary vs. continuous", " (autocorrelation = ", r, ")", sep = ""))
plot_NB = PlotPair(datapair = cbind(corrrep, Rrep[ , 8]),
         namepair = c("True latent correlation", "Estimated latent correlation"),
         title = paste("ternary vs. binary", " (autocorrelation = ", r, ")", sep = ""))
plot_NT = PlotPair(datapair = cbind(corrrep, Rrep[ , 9]),
         namepair = c("True latent correlation", "Estimated latent correlation"),
         title = paste("truncated vs. ternary", " (autocorrelation = ", r, ")", sep = ""))
plot_NN = PlotPair(datapair = cbind(corrrep, Rrep[ , 10]),
         namepair = c("True latent correlation", "Estimated latent correlation"),
         title = paste("ternary vs. ternary", " (autocorrelation = ", r, ")", sep = ""))
plot_all = gridExtra::marrangeGrob(grobs = list(plot_CC, plot_BC, plot_BB, plot_TC, plot_TB, plot_TT, plot_NC, plot_NB, plot_NT, plot_NN), widths = rep(1, 3), top = NULL, layout_matrix = matrix(c(1:12), 4, 3, byrow = T))
ggsave("serial_corr_sim.pdf", plot_all, width = 15, height = 20)
plot_con_ord = gridExtra::marrangeGrob(grobs = list(plot_CC, plot_BC, plot_BB, plot_NC, plot_NB, plot_NN), widths = rep(1, 3), top = NULL, layout_matrix = matrix(c(1:6), 2, 3, byrow = T))
ggsave("serial_corr_con_ord.pdf", plot_con_ord, width = 15, height = 10)
plot_tru = gridExtra::marrangeGrob(grobs = list(plot_TC, plot_TB, plot_TT, plot_NT), widths = rep(1, 3), top = NULL, layout_matrix = matrix(c(1:6), 2, 3, byrow = T))
ggsave("serial_corr_tru.pdf", plot_tru, width = 15, height = 10)

## BC case (AR(1))
rm(list = ls())
library(MASS)
library(mvtnorm)
library(ggplot2)
library(latentcor)
PlotPair = function(datapair, namepair = c("X", "Y"), title = "Plot X vs Y") {
  df = data.frame(datapair)
  colnames(df) = namepair
  print(ggplot(df, aes(x = datapair[ , 1], y = datapair[ , 2]))
        + geom_point(color = "blue") + geom_abline(intercept = 0, slope = 1, color = "red")
        +ggtitle(title) + xlab(namepair[1]) + ylab(namepair[2]) + xlim(-1, 1) + ylim(-1, 1))
}
corrrep = seq(-1, 1, by = .01)
Rrep = matrix(NA, length(corrrep), 2)
n = 100
rho = matrix(NA, n, n)
rho = abs(col(rho) - row(rho))
r = 0.99
rho = r^rho
p = 2
for (i in 1:length(corrrep)) {
  Z = matrix(mvrnorm(mu = rep(0, n * p), Sigma = kronecker(rho, matrix(c(1, corrrep[i], corrrep[i], 1), p, p))), ncol = p, byrow = TRUE)
# Z_1 = Z[ , 1]; Z[ , 1][Z_1 < 0] = 0; Z[ , 1][Z_1 > 0] = 1
  Z_1 = Z[ , 1]; Z_bar_1 = mean(Z_1); Z[ , 1][Z_1 <= Z_bar_1] = 0; Z[ , 1][Z_1 > Z_bar_1] = 1
  out = latentcor(Z, types = c("bin", "con"))
  tau_bin_con = out$K
  R_bin_con = out$Rpointwise[2, 1]
  re = acf(Z[,2])$acf[2]
  rhoe = matrix(NA, n, n)
  rhoe = abs(col(rhoe) - row(rhoe))
  rhoe = re^rhoe
  r_bar = (sum(rhoe) - n) / (n * (n - 1))
  zratio1 = sum(Z[ , 1] == 0) / n
  delta_1 = qnorm(zratio1)
  ff = function(s) {(4 * fMultivar::pnorm2d(delta_1, 0, rho = s/sqrt(2)) - 2 * zratio1 - sqrt(2) * s * fMultivar::dnorm2d(delta_1, 0, rho = s/sqrt(2)) * r_bar - tau_bin_con[2,1])^2}
  aug_R_bin_con = optimize(ff, lower = -0.999, upper = 0.999, tol = 1e-5)$minimum
  Rrep[i, 1] = R_bin_con; Rrep[i, 2] = aug_R_bin_con
}
PlotPair(datapair = cbind(corrrep, Rrep[ , 1]),
         namepair = c("True latent correlation", "Estimated latent correlation"),
         title = "autocorr = 0.9")
PlotPair(datapair = cbind(corrrep, Rrep[ , 2]),
         namepair = c("True latent correlation", "Augmented estimated latent correlation"),
         title = "autocorr = 0.9")

## BC case (equi-corr)
rm(list = ls())
library(MASS)
library(mvtnorm)
library(ggplot2)
library(latentcor)
PlotPair = function(datapair, namepair = c("X", "Y"), title = "Plot X vs Y") {
  df = data.frame(datapair)
  colnames(df) = namepair
  print(ggplot(df, aes(x = datapair[ , 1], y = datapair[ , 2]))
        + geom_point(color = "blue") + geom_abline(intercept = 0, slope = 1, color = "red")
        +ggtitle(title) + xlab(namepair[1]) + ylab(namepair[2]) + xlim(-1, 1) + ylim(-1, 1))
}
corrrep = seq(-.99, .99, by = .01)
Rrep = matrix(NA, length(corrrep), 2)
n = 100
rho = matrix(.7, n, n)
diag(rho) = 1
p = 2
for (i in 1:length(corrrep)) {
  Z = matrix(mvrnorm(mu = rep(0, n * p), Sigma = kronecker(rho, matrix(c(1, corrrep[i], corrrep[i], 1), p, p))), ncol = p, byrow = TRUE)
  Z_1 = Z[ , 1]; Z_bar_1 = mean(Z_1); Z[ , 1][Z_1 <= Z_bar_1] = 0; Z[ , 1][Z_1 > Z_bar_1] = 1
  out = latentcor(Z, types = c("bin", "con"))
  tau_bin_con = out$K
  R_bin_con = out$Rpointwise[2, 1]
  r_bar = acf(Z[,2])$acf[2]
  zratio1 = sum(Z[ , 1] == 0) / n
  delta_1 = qnorm(zratio1)
  ff = function(s) {(4 * fMultivar::pnorm2d(delta_1, 0, rho = s/sqrt(2)) - 2 * zratio1 - sqrt(2) * s * fMultivar::dnorm2d(delta_1, 0, rho = s/sqrt(2)) * r_bar - tau_bin_con[2,1])^2}
  aug_R_bin_con = optimize(ff, lower = -0.999, upper = 0.999, tol = 1e-5)$minimum
  Rrep[i, 1] = R_bin_con; Rrep[i, 2] = aug_R_bin_con
}
PlotPair(datapair = cbind(corrrep, Rrep[ , 1]),
         namepair = c("True latent correlation", "Estimated latent correlation"),
         title = "autocorr = 0.7")
PlotPair(datapair = cbind(corrrep, Rrep[ , 2]),
         namepair = c("True latent correlation", "Augmented estimated latent correlation"),
         title = "autocorr = 0.7")


## BB case
rm(list = ls())
library(MASS)
library(mvtnorm)
n = 1000
rho = matrix(NA, n, n)
rho = abs(col(rho) - row(rho))
r = 0.3
rho = r^rho
p = 2
Z = matrix(mvrnorm(mu = rep(0, n * p), Sigma = kronecker(rho, matrix(c(1, .5, .5, 1), p, p))), ncol = p, byrow = TRUE)
Z[ , 1][Z[ , 1] < 0] = 0
Z[ , 1][Z[ , 1] > 0] = 1
Z[ , 2][Z[ , 2] < 0] = 0
Z[ , 2][Z[ , 2] > 0] = 1
out = latentcor(Z, types = c("bin", "bin"))
tau_bin_bin = out$K
R_bin_bin = out$Rpointwise
r_bar = (sum(rho) - n) / (n * (n - 1))
zratio1 = sum(Z[ , 1] == 0) / n
delta_1 = qnorm(zratio1)
zratio2 = sum(Z[ , 2] == 0) / n
delta_2 = qnorm(zratio2)
ff = function(s) {(2 * (fMultivar::pnorm2d(delta_1, delta_2, rho = s) - zratio1 * zratio2) - 2 * s * dnorm(delta_1) * dnorm(delta_2) - tau_bin_bin[2,1])^2}
optimize(ff, lower = -0.999, upper = 0.999, tol = 1e-5)




x <- numeric(150)
x[1] <- 0
for (i in 2:150) x[i] <- 0.7*x[i-1]+rnorm(1)

y[1] <- 0
for (i in 2:150) y[i] <- 0.7*y[i-1]+rnorm(1)

x=x[51:150]
y=y[51:150]
z=x+y



## BC case (AR(1))
rm(list = ls())
library(MASS)
library(mvtnorm)
library(ggplot2)
library(latentcor)
PlotPair = function(datapair, namepair = c("X", "Y"), title = "Plot X vs Y") {
  df = data.frame(datapair)
  colnames(df) = namepair
  print(ggplot(df, aes(x = datapair[ , 1], y = datapair[ , 2]))
        + geom_point(color = "blue") + geom_abline(intercept = 0, slope = 1, color = "red")
        +ggtitle(title) + xlab(namepair[1]) + ylab(namepair[2]) + xlim(-1, 1) + ylim(-1, 1))
}
calcrho<-function(rho,rho1,rho2) {
     rho*(1-rho1*rho2)/sqrt((1-rho1^2)*(1-rho2^2))
}
burn.in<-300
n<-1000
rhosec<-seq(-1, 1, by = .01)
Rrep = matrix(NA, length(rhosec), 2)
# acf(x)
# acf(y)
for (i in 1:length(rhosec)) {
  rho<-rhosec[i]
  rho1<-0.9
  rho2<-0.9
  q12<-calcrho(rho,rho1,rho2)
  eps<-mvrnorm(n+burn.in,mu=c(0,0),Sigma=cbind(c(1,q12),c(q12,1)))
  x<-arima.sim(list(ar=rho1),n,innov=eps[burn.in+1:n,1],start.innov=eps[1:burn.in,1])
  y<-arima.sim(list(ar=rho2),n,innov=eps[burn.in+1:n,2],start.innov=eps[1:burn.in,2])
  Z = cbind(x,y)
  # Z_1 = Z[ , 1]; Z[ , 1][Z_1 < 0] = 0; Z[ , 1][Z_1 > 0] = 1
  Z_1 = Z[ , 1]; Z_bar_1 = mean(Z_1); Z[ , 1][Z_1 <= Z_bar_1] = 0; Z[ , 1][Z_1 > Z_bar_1] = 1
  out = latentcor(Z, types = c("bin", "con"))
  tau_bin_con = out$K
  R_bin_con = out$Rpointwise[2, 1]
  re = acf(Z[,2])$acf[2]
  rhoe = matrix(NA, n, n)
  rhoe = abs(col(rhoe) - row(rhoe))
  rhoe = re^rhoe
  r_bar = (sum(rhoe) - n) / (n * (n - 1))
  zratio1 = sum(Z[ , 1] == 0) / n
  delta_1 = qnorm(zratio1)
  ff = function(s) {(4 * fMultivar::pnorm2d(delta_1, 0, rho = s/sqrt(2)) - 2 * zratio1 - sqrt(2) * s * fMultivar::dnorm2d(delta_1, 0, rho = s/sqrt(2)) * r_bar - tau_bin_con[2,1])^2}
  aug_R_bin_con = optimize(ff, lower = -0.999, upper = 0.999, tol = 1e-5)$minimum
  Rrep[i, 1] = R_bin_con; Rrep[i, 2] = aug_R_bin_con
}
PlotPair(datapair = cbind(rhosec, Rrep[ , 1]),
         namepair = c("True latent correlation", "Estimated latent correlation"),
         title = "autocorr = 0.9")
PlotPair(datapair = cbind(rhosec, Rrep[ , 2]),
         namepair = c("True latent correlation", "Augmented estimated latent correlation"),
         title = "autocorr = 0.9")

ff = function(sigma) {
zratio1 = .7
de1 = qnorm(zratio1)
n = 50
rpower = 0:(n-1)
r = .5
r_all = r^rpower
sum_all = 0
for (i in 1:n) {
  sum_all = sum_all + 4 * fMultivar::pnorm2d(de1, 0, rho = (sqrt((1 - r_all[i]) / 2) * sigma))
}
return(((2 / (n * (n - 1)) * sum_all - 2 * zratio1) - .4)^2)
}
optimize(ff, lower = -0.999, upper = 0.999, tol = 1e-5)
