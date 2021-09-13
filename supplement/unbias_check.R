library(latentcor)
library(ggplot2)

rm(list=ls())
PlotPair = function(datapair, namepair = c("X", "Y"), title = "Plot X vs Y") {
  df = data.frame(datapair)
  colnames(df) = namepair
  print(ggplot(df, aes(x = datapair[ , 1], y = datapair[ , 2]))
        + geom_point(color = "blue") + geom_abline(intercept = 0, slope = 1, color = "red")
        +ggtitle(title) + xlab(namepair[1]) + ylab(namepair[2]) + xlim(-1, 1) + ylim(-1, 1))
}


# Data generation
rhos = seq(-1, 1, by = .01); rhorep = rep(NA, length(rhos)); Rrep = matrix(NA, length(rhos), 3)
for (r in 1:length(rhos)) {
  X = GenData(n = 1000, types = c("ter", "con"), rhos = rhos[r], XP = list(c(.3, .5), NA))$X
  R_nc_org = estR(X, types = c("ter", "con"), method = "original")$R
  R_nc_approx = estR(X, types = c("ter", "con"), method = "approx")$R
  rhorep[r] = rhos[r]; Rrep[r, 1] = R_nc_org[2, 1]; Rrep[r, 2] = R_nc_approx[2, 1]; Rrep[r, 3] = cor(X)[2, 1]
}

# Plot ternary/continuous case estimation via original method.
R_nc_org = PlotPair(datapair = cbind(rhorep, Rrep[,1]), namepair = c("True latent correlation", "Estimated latent correlation (original)"),
                    title = "Ternary vs. continuous")

# Plot ternary/continuous case estimation via approximation method.
R_nc_approx = PlotPair(datapair = cbind(rhorep, Rrep[,2]), namepair = c("True latent correlation", "Estimated latent correlation (approx)"),
                       title = "Ternary vs. continuous")

R_nc_pearson = PlotPair(datapair = cbind(rhorep, Rrep[,3]), namepair = c("True latent correlation", "Pearson correlation"),
                        title = "Ternary vs. continuous")
#R_nc = gridExtra::marrangeGrob(grobs = list(R_nc_org, R_nc_approx, R_nc_pearson), widths = c(10, 2, 10, 2, 10), top = NULL, layout_matrix = matrix(c(1, NA, 2, NA, 3), 1, 5))
pdf(file = "nc_org.pdf", width = 5, height = 5)
R_nc_org
dev.off()
pdf(file = "nc_approx.pdf", width = 5, height = 5)
R_nc_approx
dev.off()
pdf(file = "nc_pearson.pdf", width = 5, height = 5)
R_nc_pearson
dev.off()

cp1 = "cube"; cp2 = "cube"
for (tp1 in c("con", "bin", "ter", "tru")) {
  for (tp2 in c("con", "bin", "ter", "tru")) {
    rhorep = rep(NA, 100); Rrep = matrix(NA, 100, 4)
    for (rep in 1:100) {
      rho = runif(1, -1, 1)
      X = GenData(n = 1000, rhos = rho, types = c(tp1, tp2), copulas = c(cp1, cp2))$X
      R_org = estR(X, types = c(tp1, tp2), method = "original")$R
      R_approx = estR(X, types = c(tp1, tp2), method = "approx")$R
      R_conapprox = estR(X, types = c(ifelse(tp1 == "ter", "con", tp1), ifelse(tp2 == "ter", "con", tp2)), method = "original")$R
      rhorep[rep] = rho; Rrep[rep, 1] = R_org[2, 1]; Rrep[rep, 2] = R_approx[2, 1]; Rrep[rep, 3] = cor(X)[2, 1]; Rrep[rep, 4] = R_conapprox[2, 1]
    }
    assign(paste("R", cp1, cp2, tp1, tp2, "org", sep = "_"), Rrep[ , 1])
    assign(paste("R", cp1, cp2, tp1, tp2, "approx", sep = "_"), Rrep[ , 2])
    assign(paste("R", cp1, cp2, tp1, tp2, "pearson", sep = "_"), Rrep[ , 3])
    assign(paste("R", cp1, cp2, tp1, tp2, "conapprox", sep = "_"), Rrep[ , 4])
    PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "org", sep = "_")))),
             namepair = c("rho", paste("R", cp1, cp2, tp1, tp2, "org", sep = "_")),
             title = paste(tp1, "vs.", tp2, "(org)", sep = " "))
    PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "approx", sep = "_")))),
             namepair = c("rho", paste("R", cp1, cp2, tp1, tp2, "approx", sep = "_")),
             title = paste(tp1, "vs.", tp2, "(approx)", sep = " "))
    PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "pearson", sep = "_")))),
             namepair = c("rho", paste("R", cp1, cp2, tp1, tp2, "pearson", sep = "_")),
             title = paste(tp1, "vs.", tp2, "(pearson)", sep = " "))
    PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "conapprox", sep = "_")))),
             namepair = c("rho", paste("R", cp1, cp2, tp1, tp2, "conapprox", sep = "_")),
             title = paste(tp1, "vs.", tp2, "(conapprox)", sep = " "))
  }
}
#
# cp1 = "cube"; cp2 = "cube"
# for (tp1 in c("con")) {
#   for (tp2 in c("con", "bin", "tru", "ter", "qua", "qui", "sen", "sep", "oct", "nov", "den", "dtr")) {
#     rhorep = rep(NA, 100); Rrep = matrix(NA, 100, 3)
#     for (rep in 1:100) {
#       rho = runif(1, -1, 1)
#       X = GenData(n = 1000, rhos = rho, types = c(tp1, tp2), copulas = c(cp1, cp2))$X
#       R_org = estR(X, types = c(tp1, tp2), method = "original")$R
#       R_conapprox = estR(X, types = c("con", "con"), method = "original")$R
#       rhorep[rep] = rho; Rrep[rep, 1] = R_org[2, 1]; Rrep[rep, 2] = R_conapprox[2, 1]; Rrep[rep, 3] =  cor(X)[2, 1]
#     }
#     assign(paste("R", cp1, cp2, tp1, tp2, "org", sep = "_"), Rrep[ , 1])
#     assign(paste("R", cp1, cp2, tp1, tp2, "conapprox", sep = "_"), Rrep[ , 2])
#     assign(paste("R", cp1, cp2, tp1, tp2, "pearson", sep = "_"), Rrep[ , 3])
#     PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "org", sep = "_")))),
#              namepair = c("rho", paste("R", cp1, cp2, tp1, tp2, "org", sep = "_")),
#              title = paste(tp1, "vs.", tp2, "(org)", sep = " "))
#     PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "conapprox", sep = "_")))),
#              namepair = c("rho", paste("R", cp1, cp2, tp1, tp2, "conapprox", sep = "_")),
#              title = paste(tp1, "vs.", tp2, "(conapprox)", sep = " "))
#     PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "pearson", sep = "_")))),
#              namepair = c("rho", paste("R", cp1, cp2, tp1, tp2, "pearson", sep = "_")),
#              title = paste(tp1, "vs.", tp2, "(pearson)", sep = " "))
#   }
# }

pdf(file = "", width = , height = )
dev.off()

## Speed comparison
# p = 20
types_20 = rep(c("con", "bin", "ter", "tru"), 5)
X_20 = gen_data(types = types_20)$X
# p = 40
types_40 = rep(c("con", "bin", "ter", "tru"), 10)
X_40 = gen_data(types = types_40)$X
# p = 100
types_100 = rep(c("con", "bin", "ter", "tru"), 25)
X_100 = gen_data(types = types_100)$X
# p = 200
types_200 = rep(c("con", "bin", "ter", "tru"), 50)
X_200 = gen_data(types = types_200)$X
# p = 400
types_400 = rep(c("con", "bin", "ter", "tru"), 100)
X_400 = gen_data(types = types_400)$X

library(microbenchmark)
library(latentcor)
library(parallel)
library(doFuture)
types = list(types_20, types_40, types_100, types_200, types_400)
X = list(X_20, X_40, X_100, X_200, X_400)
md = c("original", "approx")
registerDoFuture()
plan(multicore, workers = detectCores())
timing =
foreach (i = 1:length(types), .combine = rbind) %dopar%  {
  value = rep(NA, 2)
  for (j in 1:length(md)) {
    value[j] = median(microbenchmark(latentcor(X = X[[i]], types = types[[i]], method = md[j]), times = 5L)$time)
  }
  timing_vector = value
}
df = data.frame(c("log10(20)", "log10(40)", "log10(100)", "log10(200)", "log10(400)"), log10(cbind(timing)))
print(ggplot(df, aes(x = df[ , 1])) + geom_line(y = df[ , 2], color = "blue") + geom_line(y = df[ , 3], color = "red")
      + ggtitle("Speed comparison") + xlab("log10 of dimension") + ylab("log10 of computation times"))
